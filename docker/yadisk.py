import os, json
import argparse
import urllib.parse as ul

url_temp = 'https://cloud-api.yandex.net:443/v1/disk/public/resources/download?public_key='

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Yandex disk public link downloader.')
    parser.add_argument('url', metavar='url', type=str,
                    help='Yandex disk public URL link')
    parser.add_argument('-p', '--prefix', type=str,
                        default='.',
                        help='Output prefix')
    args = parser.parse_args()
    res = os.popen('wget -qO - {}{}'.format(url_temp, args.url)).read()
    json_res = json.loads(res)
    filename = ul.parse_qs(ul.urlparse(json_res['href']).query)['filename'][0]
    os.system("wget '{}' -P '{}' -O '{}'".format(json_res['href'], args.prefix, filename))
