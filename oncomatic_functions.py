# tumor type = "BRCA" / "PAAD"
# show subtypes = yes / no
# pipeline = "mutect2" / "muse" / "varscan2" / "somaticsniper"
# number of top genes = 20
# gene of interest = show all / "ATM-BRCA1-P53"
import sys
import os.path

def readinputfile(input):
  """
  Reads the input file (which must be formatted correctly). It then extracts the
  params inside the input file and executes the appropriate script.
  """
  # check if file exists - exits with error message if file cannot be found
  _ensure_file_exists(input)
  in_file  = open(input, "r")

  for words in in_file.readlines():
    print(words)
    in_file.close()
    
def _ensure_file_exists(file):
  """
  Simply check if the input file exists. If the file is not found, then shows
  an error message and then exits.
  """
  if os.path.isfile(file) == False:
    print "####################################################################"
    print "# IOError: The following input file was not found:                 #"
    print "# ", file, "                                                       #"
    print "# Please try running oncomatic.py again but ensure that all        #"
    print "# filenames are correct.                                           #"
    print "####################################################################"
    sys.exit()


def help(help_file="help.txt"):
  """
  Prints the help text from the help_file and then exits
  """
  with open(help_file) as f:
    print # a blank line
    print f.read()
    sys.exit()

def menu():
  """
  Runs the menu feature - i.e. the interactive usage
  """
  menu_text()
  option = raw_input('Please pick an option: ')
  if option == 'help':
    help()
  elif option == '1':
    print # a blank line
    first_file  = raw_input("Please enter the (full or relative) path to your input File: ")
    readinputfile(first_file)
  elif option == '2':
      help()
  elif option == '3':
    print # a blank line
    print "####################################################################"
    print "# Goodbye! Thank you for using oncomatic. Please cite responsibly. #"
    print "####################################################################"
    print # a blank line
    sys.exit()
  else:
    print # a blank line
    print "####################################################################"
    print "# You have entered an unrecognised term. Please re-read the menu   #"
    print "# options and try again. Type 'help' for a readme.                 #"
    print "####################################################################"  
    print # a blank line
    print # a blank line
    menu() # Show the menu again.

def menu_text():
  """
  Display the menu text
  """
  print # a blank line
  print "MENU"
  print "----"
  print '[1] Print a dotplot of two input FASTA sequences.'
  print '[2] or [help] Bring up help instructions.'
  print '[3] Quit the program.'
  print # a blank line