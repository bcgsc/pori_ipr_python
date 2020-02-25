import hashlib
import json


def hash_key(key):
    body = json.dumps({'key': key}, sort_keys=True)
    hash_code = hashlib.md5(body.encode('utf-8')).hexdigest()
    return hash_code
