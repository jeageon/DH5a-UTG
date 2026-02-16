from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Union
from urllib.parse import urlencode

import requests
from pydantic import BaseModel

from .exceptions import ToolError
from ..config import USER_AGENT


class ResponseWrapper(BaseModel):
    url: str
    status_code: int
    headers: Dict[str, str] = {}
    text: Optional[str] = None
    json_obj: Optional[Union[Dict[str, Any], List[Any]]] = None


class ApiClient:
    def __init__(
        self,
        timeout: float = 20.0,
        retries: int = 5,
        backoff_factor: float = 1.0,
        cache_enabled: bool = True,
        cache_path: Optional[Union[str, Path]] = None,
        ttl_hours: int = 24,
        offline: bool = False,
        user_agent: str = USER_AGENT,
    ) -> None:
        self.timeout = timeout
        self.retries = retries
        self.backoff_factor = backoff_factor
        self.cache_enabled = cache_enabled
        self.offline = offline
        self.ttl_seconds = ttl_hours * 3600
        self.cache_path = Path(cache_path or Path("data/cache"))
        self.cache_file = self.cache_path / "utg_api_cache.json"
        if self.cache_enabled:
            self.cache_path.mkdir(parents=True, exist_ok=True)
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": user_agent})

    def _build_cache_key(
        self,
        method: str,
        url: str,
        params: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None,
        data: Any = None,
        json_payload: Any = None,
    ) -> str:
        method = method.upper()
        payload = {
            "method": method,
            "url": url,
            "params": params or {},
            "headers": {**self.session.headers, **(headers or {})},
            "data": data,
            "json": json_payload,
        }
        body = json.dumps(payload, sort_keys=True, default=str)
        return hashlib.sha256(body.encode("utf-8")).hexdigest()

    def _read_cache(self) -> dict[str, Any]:
        if not self.cache_file.exists():
            return {}
        try:
            text = self.cache_file.read_text(encoding="utf-8")
            raw = json.loads(text)
            return raw if isinstance(raw, dict) else {}
        except Exception:
            return {}

    def _write_cache(self, payload: dict[str, Any]) -> None:
        if not self.cache_enabled:
            return
        try:
            self.cache_file.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        except Exception:
            pass

    def _get_cached(self, key: str) -> Optional[ResponseWrapper]:
        if not self.cache_enabled:
            return None
        data = self._read_cache()
        cached = data.get(key)
        if not cached:
            return None
        saved = cached.get("saved_at", 0)
        if self.ttl_seconds > 0 and (time.time() - saved) > self.ttl_seconds:
            return None
        return ResponseWrapper(
            url=cached.get("url", ""),
            status_code=cached.get("status_code", 0),
            headers=cached.get("headers", {}),
            text=cached.get("text"),
            json_obj=cached.get("json_obj"),
        )

    def _cache_if_needed(self, key: str, response: ResponseWrapper) -> ResponseWrapper:
        if not self.cache_enabled:
            return response
        data = self._read_cache()
        data[key] = {
            "url": response.url,
            "status_code": response.status_code,
            "headers": response.headers,
            "text": response.text,
            "json_obj": response.json_obj,
            "saved_at": time.time(),
        }
        self._write_cache(data)
        return response

    def _parse_response(self, url: str, response: requests.Response) -> ResponseWrapper:
        text = response.text
        parsed = None
        ctype = response.headers.get("content-type", "").lower()
        if "json" in ctype or response.text.lstrip().startswith("{") or response.text.lstrip().startswith("["):
            try:
                parsed = response.json()
            except Exception:
                parsed = None
        return ResponseWrapper(
            url=url,
            status_code=response.status_code,
            headers=dict(response.headers),
            text=text,
            json_obj=parsed,
        )

    def _sleep(self, attempt: int, retry_after: Optional[float] = None) -> None:
        if attempt <= 0:
            return
        if retry_after is not None:
            delay = max(0.0, retry_after)
        else:
            delay = self.backoff_factor * (2 ** (attempt - 1))
        delay += min(1.0, 0.25 * attempt)
        time.sleep(delay)

    def _request(
        self,
        method: str,
        url: str,
        headers: Optional[Dict[str, str]] = None,
        params: Optional[Dict[str, Any]] = None,
        data: Any = None,
        json_payload: Any = None,
        disable_cache: bool = False,
    ) -> ResponseWrapper:
        key = self._build_cache_key(
            method,
            url,
            params=params,
            headers=headers,
            data=data,
            json_payload=json_payload,
        )
        cached = None if disable_cache else self._get_cached(key)
        if cached:
            return cached

        if self.offline:
            raise ToolError(f"Offline mode: cache miss for {method} {url} {params}")

        merged_headers = {**self.session.headers, **(headers or {})}

        for attempt in range(self.retries + 1):
            try:
                response = self.session.request(
                    method=method,
                    url=url,
                    headers=merged_headers,
                    params=params,
                    data=data,
                    json=json_payload,
                    timeout=self.timeout,
                )
            except requests.RequestException as exc:
                if attempt < self.retries:
                    self._sleep(attempt + 1)
                    continue
                raise ToolError(f"Network error for {url}: {exc}") from exc

            if response.status_code == 429:
                if attempt < self.retries:
                    retry_after = response.headers.get("Retry-After")
                    delay = float(retry_after) if retry_after else None
                    self._sleep(attempt + 1, delay)
                    continue
                raise ToolError(f"Rate limit hit for {url}: {response.status_code}")

            if 500 <= response.status_code < 600 and attempt < self.retries:
                self._sleep(attempt + 1)
                continue

            if response.status_code >= 400:
                message = response.text[:500]
                raise ToolError(
                    f"Request failed ({response.status_code}) for {url} {urlencode(params or {}, doseq=True)}: {message}"
                )

            parsed = self._parse_response(url, response)
            if disable_cache:
                return parsed
            return self._cache_if_needed(key, parsed)

        raise ToolError(f"Request exhausted retries for {url}")

    def get(
        self,
        url: str,
        headers: Optional[Dict[str, str]] = None,
        params: Optional[Dict[str, Any]] = None,
        disable_cache: bool = False,
    ) -> ResponseWrapper:
        return self._request("GET", url, headers=headers, params=params, disable_cache=disable_cache)

    def post(
        self,
        url: str,
        headers: Optional[Dict[str, str]] = None,
        data: Any = None,
        json_payload: Any = None,
        disable_cache: bool = False,
    ) -> ResponseWrapper:
        return self._request(
            "POST",
            url,
            headers=headers,
            data=data,
            json_payload=json_payload,
            disable_cache=disable_cache,
        )
