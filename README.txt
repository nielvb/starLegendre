Code for star-framework: solving non-autonomous linear systems of ODEs


--------------- Files ---------------------
- 'Init.m' -> loads in code necessary to run the algorithms, note that external libraries (chebfun) is required, for more information see 'Init.m'.

- 'genCoeffMatrix.m' -> computes a matrix which contains the Legendre coefficients of a smooth function f(t) multiplied by the (bivariate) Heaviside step function theta(t-s)

- 'StarLegendre_scalar.m' -> Solves the scaler ODE u'(t) = f(t) u(t), u(-1)=1, on [-1,1] via the *-approach.

- 'StarLegendre_matrix.m' -> Solves the matrix ODE U'(t) = H(t) U(t), U(-1)=I, on [-1,1] via the *-approach.


---------------- Folders -------------------------
- 'Experiments' -> Contains illustrative examples how to use the software.

- 'Papers' -> Contains the numerical experiments performed in the papers (see below for a list of papers).

- Note: each folder has a dedicated Init.m to load in additional code.

---------------- Papers ---------------------------
[1] https://arxiv.org/pdf/2303.11284.pdf (to be published in ETNA)

[2] https://arxiv.org/abs/2311.04144



---------------- Copyright Statement-----------------------

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: 

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. 

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
