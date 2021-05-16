// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use dusk_bls12_381::BlsScalar as Scalar;

// Apply affine transformation to state via MDS matrix multiplication
pub fn concrete(
    state: [Scalar; 3],
    matrix: [[Scalar; 3]; 3],
    constants: [Scalar; 3],
) -> [Scalar; 3] {
    let mut new_state = constants;

    // matrix multiplication
    for i in 0..3 {
        for j in 0..3 {
            new_state[i] += matrix[i][j] * state[j];
        }
    }

    new_state
}

mod tests {
    use super::*;
    use crate::constants::{CONSTANTS_BLS, MATRIX_BLS};

    #[test]
    fn test_concrete() {
        let state = [Scalar::from(4), Scalar::from(3), Scalar::from(2)];
        let matrix = MATRIX_BLS;
        let constants = CONSTANTS_BLS;
        let output = concrete(state, matrix, constants[0]);

        let copy_matrix = [
            [Scalar::from(2), Scalar::from(1), Scalar::from(1)],
            [Scalar::from(1), Scalar::from(2), Scalar::from(1)],
            [Scalar::from(1), Scalar::from(1), Scalar::from(2)],
        ];

        let mut new_state = constants[0];
        for i in 0..3 {
            for j in 0..3 {
                new_state[i] += copy_matrix[i][j] * state[j];
            }
        }

        assert_eq!(new_state, output);
    }
}
