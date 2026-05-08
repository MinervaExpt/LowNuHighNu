#!/usr/bin/env python3

"""Renew the token for accessing the MINERvA vault"""

# H. Schellman 9-20-2025
# This does a token renewal, if you have not gotten a vault token yet, it prompts you to get one.
# argument is an optional override of the 600 second sleep time.
# vscode copilot made suggestions which I did not always follow.

import sys, subprocess, time, datetime, socket

def renew_token(log):
    "" 
    cmd = "htgettoken -q -i minerva --vaultserver htvaultprod.fnal.gov --nooidc".split(" ")
    now = datetime.datetime.now().strftime("%Y-%m-%d-%H%M.%S")
    try:
        log.write ("%s Renewing token with command: %s\n"%(now," ".join(cmd)))
        subprocess.run(cmd, check=True)
    
        log.write("%s Token renewed\n"%now)
    except:
        log.write("%s WARNING: Token renewal failed, skip for now\n"%now)
        print ("WARNING: Token renewal failed: You might want to try doing:\nhtgettoken -i minerva --vaultserver htvaultprod.fnal.gov \nto get a fresh vault token")
    log.flush()

if __name__ == "__main__":
    alive = 0
    step = 600
    # override step if given as argument
    # also see if someone wanted help
    log = open("renew_token_%s.log"%socket.gethostname(),"w")
    if len(sys.argv)>1:
        if "h" in sys.argv[1].lower():
            print ("Usage: %s [step_seconds]"%sys.argv[0])
            print ("This will renew your MINERvA vault token every step_seconds (default 600 seconds)")
            sys.exit(0)
        step = sys.argv[1]
        step = int(step)
    while alive < 3600*24: # run for 12 hours
        renew_token(log)
        time.sleep(step) # sleep for 600 seconds
        alive += step
    log.close()
