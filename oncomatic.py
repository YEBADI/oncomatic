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
    elif len(sys.argv) == 3:
      # if only one sys.arv after the python script name
      print # a blank line
      print "##################################################################"
      print "# We noticed that you entered more than one file, oncoprint only #"
      print "# handles one input file.                                        #"
      print "##################################################################"
      print # a blank line
    else:
      # continue with analysis
      first_file  = sys.argv[1]
      readinputfile(first_file)
  else:
    menu() # start the menu function instead.

# Start the Applet
run_oncomatic()