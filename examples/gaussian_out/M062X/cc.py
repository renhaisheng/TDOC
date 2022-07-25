#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from re import findall

with open('[CH2]-T.out') as f:
	content=f.read()
res = list(map(float, findall('before annihilation\s+(\S+?),\s+after\s+(\S+)\n', content)[-1]))
print(res)