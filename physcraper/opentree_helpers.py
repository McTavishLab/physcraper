import requests
import json
import sys

if sys.version_info < (3,):
    from urllib2 import HTTPError
#    from urllib2 import ConnectionError
else:
    from urllib.error import HTTPError


def check_if_ottid_in_synth(ottid):
    url = 'https://api.opentreeoflife.org/v3/tree_of_life/node_info'
    payload = json.dumps({"ott_id":int(ottid)})
    headers = {'content-type': 'application/json', 'Accept-Charset': 'UTF-8'}
    try:
        r = requests.post(url, data=payload, headers=headers)
        if r.status_code == 200:
            return 1
        elif r.status_code == 400:
            return 0
        elif r.status_code == 502:
            sys.stderr.write("Bad OpenTree taxon ID: {}".format(ottid))
            return 0
        else:
            sys.stderr.write("unexpected status code from node_info call: {}".format(r.status_code))
            return 0
    except requests.ConnectionError:
        sys.stderr.write("Connection Error - coud not get taxon information from OpenTree\n")
