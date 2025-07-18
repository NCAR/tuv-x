#!/usr/bin/env python3
# utilsLite.py
#
# Minimal utilities for Python scripts.
#
# Created by Carl Drews - July 2020

import sys
import re



# Sanitize input coming from the browser.
def sanitize(original):
	if (original == None):
		return ""

	# strip everything but alphanumeric chars, colon, underscore, and hyphen
	safeString = re.sub(r'[^a-zA-Z0-9.:_-]+', '', original)
	return safeString



# Convert safely from string to integer (alphas convert to 0).
def safeInt(intString):
	intValue = 0

	try:
		intValue = int(intString)
	except ValueError as error:
		intValue = 0

	return intValue



# Convert string to number, or 0.0 if not numeric.
# numString = string that probably can be converted
def safeFloat(numString):
   result = -1.0
   try:
      result = float(numString)
   except ValueError:
      result = 0.0

   return result



# pass PNG image back to the browser
# response = django.http.HttpResponse object containing binary image
def imageToBrowser(response):
   # send the http header
   for item in response.items():
      print("{}: {}".format(item[0], item[1]))
   print()
   sys.stdout.flush()

   # send the binary contents of the image
   sys.stdout.buffer.write(response.content)
   sys.stdout.flush()

   return

