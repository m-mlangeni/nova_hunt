"""
tns_update.py -- download the current TNS public-objects list with the MeerLICHT/BlackGEM TNS bot.

  python tns_update.py                                  # key from api_key_MeerLICHT_Bot1_TNS.txt
  python tns_update.py --key-file /path/to/key.txt
  TNS_API_KEY=... python tns_update.py                  # or from an environment variable

Writes ../BG_Trans_Query_ID/tns_public_objects.csv (the file nova_hunt.py reads by default),
keeping the previous version as tns_public_objects.csv.bak. The first line of the TNS file
(the object count) is dropped so the CSV starts with the column names.

The API key is a credential: keep the key file out of git (it is in .gitignore) and out of chats.
"""
import argparse
import io
import os
import sys
import zipfile

import requests

URL = "https://www.wis-tns.org/system/files/tns_public_objects/tns_public_objects.csv.zip"
BOT = 'tns_marker{"tns_id":59294, "type":"bot", "name":"MeerLICHT_Bot1"}'


def read_key(key_file):
    key = os.environ.get("TNS_API_KEY", "").strip()
    if key:
        return key
    if not os.path.exists(key_file):
        sys.exit(f"No TNS API key: set TNS_API_KEY or put the key in {key_file}")
    with open(key_file) as fh:
        return fh.read().strip()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--key-file", default="api_key_MeerLICHT_Bot1_TNS.txt")
    ap.add_argument("--out", default="../BG_Trans_Query_ID/tns_public_objects.csv")
    a = ap.parse_args()

    key = read_key(a.key_file)
    print("Downloading the TNS public-objects list ...")
    r = requests.post(URL, headers={"user-agent": BOT}, data={"api_key": key}, timeout=600)
    if r.status_code != 200:
        sys.exit(f"TNS refused the download: HTTP {r.status_code} {r.text[:200]}")
    try:
        z = zipfile.ZipFile(io.BytesIO(r.content))
    except zipfile.BadZipFile:
        sys.exit(f"TNS did not return a zip file (first bytes: {r.content[:200]!r}); "
                 "check the API key and bot name")
    name = next(n for n in z.namelist() if n.endswith(".csv"))
    lines = z.read(name).decode("utf-8").splitlines(keepends=True)
    if not lines or "objid" not in lines[1]:
        sys.exit("Unexpected TNS file layout: second line should hold the column names")

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    if os.path.exists(a.out):
        os.replace(a.out, a.out + ".bak")
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.writelines(lines[1:])
    print(f"Saved {len(lines) - 2} objects to {a.out} (TNS header line: {lines[0].strip()})")


if __name__ == "__main__":
    main()
