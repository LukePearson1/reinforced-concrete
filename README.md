# Reinforced Concrete
Implementation of the Hash function, Reinforced 
Concrete, for BLS12-381 curve group. This implementation 
is designed for zero-knowledge in-circuit hashing and makes use of 
Plookup tables to lookup values,
Implementation designed by the 
[dusk](https://dusk.network) team. 

## About
This hash function was developed by Dmitry 
Khovratovich et al. and makes use of lookup 
tables to greatly reduce the constraint count as 
compared to previous state of the art zero knowledge hash functions.
The website showing the asscoiated 
[paper](https://drive.google.com/file/d/1MCIqD8XwKrHVBQPc26XjAmM9RyrLDpjw/view) and reference 
implementation can be found [here](https://www.rc-hash.info/). 

RC takes in a state of three elements from the 
desired field F_p, and outputs the hash, three 
elements from the same F_p. In practice RC will 
be used in a sponge construction, meaning one 
can input an arbitrary number of field elements 
and get back as many elements as desired. Additional 
information about this sponge construction can be 
found [here](https://en.wikipedia.org/wiki/Sponge_function).

The security of this hash function is also not 
based solely on algebraic degree, and uses more 
traditional components that have received decades
of scrutiny - the security of this hash function 
is thus regarded as higher than many of its 
competitors.

**This library is as follows:**


RC is composed of three functions: concrete, 
bricks and bar. Each of these functions contributes 
to different security guarantees of the hash function.
Each function will be explained more below; the figure
below illustrates the order in which these functions
are applied.

![construction](https://user-images.githubusercontent.com/49643572/129221137-b68ad83f-1cdb-4643-bb50-302aa87bc3ac.png)


*Concrete*
This function multiplies the three element state by a 
3x3 MDS matrix, and consequently adds a constant to each
of the three elements. This is displayed concisely in 
the figure below. Note that the index $j$ indicates the
round number (concrete is used six times), so the three
dimensional constant vector $c^(j)$ is different for
each of these rounds.

![concrete](https://user-images.githubusercontent.com/49643572/129221363-f852674d-3180-499b-a715-718d9d773406.png)


Concrete allows statistical and algebraic properties 
to be spread to the whole state, providing protection.

*Bricks*
Bricks is a non-linear permutation of the three entry 
state, where the permutation is decided by the smallest
$d$ that satisfies $gcd(p-1, d) = 1. Where $p$ is the size of the field.

This repo uses [BLS12_381](https://github.com/dusk-network/bls12_381/tree/master/src), so $d$ takes the value 5 in the associated repo.

*Bars*
This is the aspect of the hash function that differntiates 
the methodology from previous constructions. 
By utilising lookups, it greatly reduces the constraint 
count compared to previous functions. The function 
begins by decomposing a single element of F_p; as 
opposed to working on all three elements at once,
this function is merely applied to each entry 
separately.
 
Numbers can be decomposed according to a basis. 
For example, 4785 can be separated into 4 buckets,
(4,7,8,5), which is according to the basis 
(10^3, 10^2, 10^1, 10^0). We can actually use 
any desired basis for this. Bar uses a basis 
consisting of $n$ elements; in the case of BLS12_381,
$n=27$. (The choice of n, and the subsequent 
basis elements is too detailed for this overview).
This gives rise to the basis 
(693, 696, 694, 668, 679, 695, 691, 693, 700, 688, 700, 694, 701, 694, 699, 701, 701, 701, 695, 698, 697, 703, 702, 691, 688, 703, 679) for BLS12_381. Note that the order does 
matter. An input value $x$ will decomposed 
into 27 elements (x_1, x_2, ..., x_27) such 
that $x = x_1*s_2*s_3*...*s_27 + x_2*s_3*s_4*...*s_27 + ... + x_26*s_27 + x_27$, where s_2 = 696, s_3 = 694 
and so on. Note that it is not a mistake that 
s_1 does not appear in the composition.
We are first interested in the decomposition 
of the largest element in the field, $p-1$. 
It will be decomposed into something of the 
form (v_1, v_2, ..., v_27), where each $v_i < s_i$. The largest prime that is smaller than each v_i is 659 - we 
will use this to determine the permutation. 
This reason for this choice is too in depth 
for this overview, but more detail can be 
found in the associated paper 
(the short of it is that it ensures collision resistance and overflow).
Now when we decompose $x$ and get the 27 elements, 
we permute each $x_i$ accordingly: if $x_i \geq 659$,
then it does not change. If $x_i < 659$, 
then it is permuted according to some 
fixed function (i.e. whenever RC is used
on a particular field, this function is 
always the same, so 0 is always sent to 
15 on BLS12_381). This permutation uses 
a plookup gate and takes only one constraint.
Once each of these 27 items have been permuted,
we have a new 27-tuple (y_1, y_2, ..., y_27). We then recompose this tuple in the reverse way to how the 
decomposition was done; i.e. we compute 
$y = y_1*s_2*s_3*...*s_27 + y_2*s_3*s_4*...*s_27 + ... + y_26*s_27 + x_27$.
Alongside this are various constraints 
that ensure that everything was done in a valid manner.

## Licensing
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Copyright (c) DUSK NETWORK. All rights reserved.
