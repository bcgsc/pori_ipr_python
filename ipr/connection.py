import requests

import json
import os
import zlib
from typing import Dict, List

from .constants import DEFAULT_URL

IMAGE_MAX = 20  # cannot upload more than 20 images at a time


class IprConnection:
    def __init__(self, username: str, password: str, url: str = DEFAULT_URL):
        self.token = None
        self.url = url
        self.username = username
        self.password = password
        self.headers = {
            'Accept': 'application/json',
            'Content-Type': 'application/json',
            'Content-Encoding': 'deflate',
        }
        self.cache: Dict[str, List[Dict]] = {}
        self.request_count = 0

    def request(self, endpoint: str, method: str = 'GET', **kwargs) -> Dict:
        """Request wrapper to handle adding common headers and logging

        Args:
            endpoint (string): api endpoint, excluding the base uri
            method (str, optional): the http method. Defaults to 'GET'.

        Returns:
            dict: the json response as a python dict
        """
        url = f'{self.url}/{endpoint}'
        self.request_count += 1
        headers = kwargs.pop('headers', self.headers)
        resp = requests.request(
            method, url, headers=headers, auth=(self.username, self.password), **kwargs
        )
        try:
            resp.raise_for_status()
        except requests.exceptions.HTTPError as err:
            # try to get more error details
            message = str(err)
            try:
                message += ' ' + resp.json()['error']['message']
            except Exception:
                pass

            raise requests.exceptions.HTTPError(message)
        return resp.json()

    def post(self, uri: str, data: Dict = {}, **kwargs) -> Dict:
        """Convenience method for making post requests"""
        return self.request(
            uri,
            method='POST',
            data=zlib.compress(json.dumps(data, allow_nan=False).encode('utf-8')),
            **kwargs,
        )

    def upload_report(self, content: Dict) -> Dict:
        return self.post('/reports', content)

    def set_analyst_comments(self, report_id: str, data: Dict) -> Dict:
        """
        Update report comments to an existing report

        TODO:
            Add to main upload.
            Pending: https://www.bcgsc.ca/jira/browse/DEVSU-1177
        """
        return self.request(
            f'/reports/{report_id}/summary/analyst-comments',
            method='PUT',
            data=zlib.compress(json.dumps(data, allow_nan=False).encode('utf-8')),
        )

    def post_images(self, report_id: str, files: Dict[str, str], data: Dict[str, str] = {}) -> None:
        """
        Post images to the report
        """
        file_keys = list(files.keys())
        start_index = 0
        image_errors = set()
        while start_index < len(file_keys):
            current_files = {}
            for key in file_keys[start_index : start_index + IMAGE_MAX]:
                path = files[key]
                if not os.path.exists(path):
                    raise FileNotFoundError(path)
                current_files[key] = path
            open_files = {k: open(f, 'rb') for (k, f) in current_files.items()}
            try:
                resp = self.request(
                    f'reports/{report_id}/image',
                    method='POST',
                    data=data,
                    files=open_files,
                    headers={},
                )
                for status in resp:
                    if status.get('upload') != 'successful':
                        image_errors.add(status['key'])
            finally:
                for handler in open_files.values():
                    handler.close()
            start_index += IMAGE_MAX
        if image_errors:
            raise ValueError(f'Error uploading images ({", ".join(sorted(list(image_errors)))})')
