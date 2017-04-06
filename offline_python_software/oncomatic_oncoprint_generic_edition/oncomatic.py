# import all functions from functions file
from oncomatic_functions import *

def run_oncomatic():
  """
  User can run oncomatic and get menu and do file input or can overide menu and
  simply enter "python oncomatic.py [FILE 1]".
  This code checks if the user has decided to place the file in command line
  or if they have decided to simply initiate the menu.
  """

  if len(sys.argv) > 1: 
    if 'help' in sys.argv:
      # if the string 'help is shown anywhere in the system arguments, then '
      # show the help text and exit.
      help()
      sys.exit()
    elif len(sys.argv) >= 7:
      # if only one sys.arv after the python script name
      print # a blank line
      print "##################################################################"
      print "# We noticed that you entered more than 5 terms, oncoprint only  #"
      print "# handles 5 arguments, please write 'python oncomatic.py help'   #"
      print "# for a readme file to understand how to use oncomatic.          #"
      print "##################################################################"
      print # a blank line
    else:
      # continue with analysis
      arg1  = sys.argv[1] #tumor type
      arg2  = sys.argv[2] #pipeline
      arg3  = sys.argv[3] #number of genes to show
      arg4  = sys.argv[4] #show specific genes?
      arg5  = sys.argv[5] #if yes, state them
      makeonco(arg1, arg2, arg3, arg4, arg5)

  else:
    menu() # start the menu function instead.

# Start the Applet
run_oncomatic()