# /******************************************************************************
# * Copyright (c) 2023 Thales
# *
# * Permission is hereby granted, free of charge, to any person obtaining a copy
# * of this software and associated documentation files (the "Software"), to deal
# * in the Software without restriction, including without limitation the rights
# * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# * copies of the Software, and to permit persons to whom the Software is
# * furnished to do so, subject to the following conditions:
# * *
# * The above copyright notice and this permission notice shall be included in
# * all copies or substantial portions of the Software.
# * *
# * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# * THE SOFTWARE.
# -
# ******************************************************************************/




SHELL := /bin/bash # Use bash syntax
CC := g++
LDFLAGS1 := -mpclmul
LDFLAGS2 := -Lfichiers_complementaires/keccak/build -lkeccak
output := helium-rmfe_ascon

all:
	@${CC} -std=c++17 test_helium-rmfe_ascon.cpp -o ${output} ${LDFLAGS1} ${LDFLAGS2}

	@# number_rep : nombre de répétition que l'on fait du test.
	@ number_rep=1 ; number=1 ; while [[ $$number -le number_rep ]] ; do \
		./helium-rmfe_ascon; \
		((number = number + 1)) ; \
		# sleep 3 ; \
	done
