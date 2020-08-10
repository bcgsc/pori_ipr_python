import json
import zlib
from typing import Dict, List

import requests

from .constants import DEFAULT_URL


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
        resp = requests.request(
            method, url, headers=self.headers, auth=(self.username, self.password), **kwargs
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
