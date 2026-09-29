"""Direct LAN-only XINAN transport. No global config changes or proxy."""
import base64,subprocess,pathlib,json,argparse
HOST='USER@LAN_HOST'
KEY=str(pathlib.Path.home()/'.ssh/A15_IDENTITY_FILE')
REMOTE='D:/AI/A15_BRIGHT_DARK_FULLWAVE'
OPTIONS=['-i',KEY,'-o','BatchMode=yes','-o','ConnectTimeout=8','-o','ProxyCommand=none','-o','ProxyJump=none']
def powershell(script,timeout=60):
    script='$ProgressPreference="SilentlyContinue"; [Console]::OutputEncoding=[System.Text.Encoding]::UTF8; '+script
    encoded=base64.b64encode(script.encode('utf-16le')).decode()
    r=subprocess.run(['ssh',*OPTIONS,HOST,'powershell -NoProfile -EncodedCommand '+encoded],capture_output=True,text=True,timeout=timeout)
    if r.returncode:raise RuntimeError(r.stdout+'\n'+r.stderr)
    return r.stdout

def put(local,remote):
    subprocess.run(['scp',*OPTIONS,str(local),HOST+':'+remote],check=True,capture_output=True)

def get(remote,local):
    subprocess.run(['scp',*OPTIONS,HOST+':'+remote,str(local)],check=True,capture_output=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--script');p.add_argument('--timeout',type=int,default=60);a=p.parse_args();print(powershell(a.script,a.timeout))
