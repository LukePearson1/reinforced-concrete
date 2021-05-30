// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use crate::constants::{MONTGOMERY_FOUR, MONTGOMERY_THREE, MONTGOMERY_TWO};
use dusk_plonk::bls12_381::BlsScalar as Scalar;

// Element-wise power function
// α1 = 1
// α2 = 3
// β1 = 2
// β2 = 4
// The above constants are used in the
// paper as part of the bricks description.
// These are given in a montgomery reduction
// form
pub fn brick(state: [Scalar; 3]) -> [Scalar; 3] {
    let mut new_state = [Scalar::zero(); 3];

    let x_squared = state[0] * state[0];
    // From the description of alpha_i - (4 * beta_i) != a square modulo p
    // d is taken to be 5. Thus x1^5 is the first element in state output.
    new_state[0] = x_squared * x_squared * state[0];
    new_state[1] = state[1] * (&x_squared + state[0] + &MONTGOMERY_TWO);
    new_state[2] = state[2]
        * ((state[1] * state[1])
            + MONTGOMERY_THREE * state[1]
            + MONTGOMERY_FOUR);
    new_state
}

mod tests {
    use super::*;

    #[test]
    fn test_bricks() {
        let input = [Scalar::from(4), Scalar::from(3), Scalar::from(2)];
        let output = brick(input);

        let two = Scalar([
            17179869180,
            12756850513266774020,
            3681868479150465002,
            3479420709561305823,
        ]);
        let three = Scalar([
            25769803770,
            688531696190609414,
            14746174755580473312,
            5219131064341958734,
        ]);
        let four = Scalar([
            34359738360,
            7066956952823996424,
            7363736958300930005,
            6958841419122611646,
        ]);
        let quadratic_x = Scalar::from(4) * Scalar::from(4);
        let quadratic_y = Scalar::from(3) * Scalar::from(3);
        let element_1 = quadratic_x * quadratic_x * Scalar::from(4);
        let element_2 = Scalar::from(3) * (quadratic_x + Scalar::from(4) + two);
        let element_3 =
            Scalar::from(2) * (quadratic_y + (three * Scalar::from(3)) + four);
        let calculated_output = [element_1, element_2, element_3];
        assert_eq!(output[0], calculated_output[0]);
        assert_eq!(output[1], calculated_output[1]);
        assert_eq!(output[2], calculated_output[2]);
    }
}
