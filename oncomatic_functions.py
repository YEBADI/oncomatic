# tumor type = "BRCA" / "PAAD"
# show subtypes = yes / no
# pipeline = "mutect2" / "muse" / "varscan2" / "somaticsniper"
# number of top genes = 20
# gene of interest = show all / "ATM-BRCA1-P53"
import sys
import os.path
import subprocess

def makeonco(arg1, arg2, arg3, arg4, arg5):
  subprocess.call(['Rscript', 'all_tumors_oncoprint', arg1, arg2, arg3, arg4, arg5])


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
    arg1  = raw_input("Please enter TCGA tumor type (e.g. 'BRCA' or 'PAAD' or 'LUAD'): ")
    arg2  = raw_input("Please write pipeline: 'mutect2', 'varscan2', 'muse', or 'somaticsniper': ")
    arg3  = raw_input("Please enter number of genes to show (e.g. '20' or '40'): ")
    arg4  = raw_input("Would you like to only show specific genes? Write 'yes' or 'no': ")
    arg5  = raw_input("If 'no', press RETURN KEY to skip. If 'yes' please write out genes in capitals seperated by a dash (e.g. 'ATM-P53-BRCA1'): ")

    makeonco(arg1, arg2, arg3, arg4, arg5)
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
  print '[1] Generate an oncoprint for top genes or specific genes of choice.'
  print '[2] or [help] Bring up help instructions.'
  print '[3] Quit the program.'
  print # a blank line