#!/usr/bin/env python

import os
import sys

from parse import *

input = sys.stdin.readline()
result = parse("RESULT: {res} (minimum: {min}%, actual: {act}%)", input.strip())
minimum = float(result.named["act"]) - float(os.environ["THRESHOLD"])

print(round(minimum,2))

sys.exit(0)
