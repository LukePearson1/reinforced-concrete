// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! This module contains the fundamental functions required for hashing,
//! using lookups. These are the three main functions of Zelbet:
//! Bricks, Concrete and Bars.

extern crate dusk_bls12_381 as BLS;

mod bar;
mod concrete;
mod brick;

use BLS::BlsScalar as Scalar;
pub use concrete::concrete;
pub use bar::bar;
pub use brick::brick;

/// Reinforced concrete hash function, taking in the hash parameters and
/// three-element item to be hashed, and outputting the hash value (three BLS
/// scalar elements)
pub fn zelbet_hash(
    scalar_inputs: [Scalar; 3],
    matrix: [[Scalar; 3]; 3],
    constants: [[Scalar; 3]; 6],
) -> [Scalar; 3] {
    let mut new_state = concrete(scalar_inputs, matrix, constants[0].clone());
    new_state = brick(new_state);
    new_state = concrete(new_state, matrix, constants[1].clone());
    new_state = brick(new_state);
    new_state = concrete(new_state, matrix, constants[2].clone());
    bar(&mut new_state);
    new_state = concrete(new_state, matrix, constants[3].clone());
    new_state = brick(new_state);
    new_state = concrete(new_state, matrix, constants[4].clone());
    new_state = brick(new_state);
    new_state = concrete(new_state, matrix, constants[5].clone());
    new_state
}
