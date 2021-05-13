// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

// This file contains the circuit implementation of the
// zelbet hash function

use dusk_plonk::constraint_system::Variable;
use dusk_plonk::prelude::*;

/// xd1,x2 ·(x21 +α1 ·x1 +β1),x3 ·(x2 +α2 ·x2 +β2))

/// This function computes the in-circuit brick function,
/// as part of the hashing gadget
pub fn brick_gagdet(
    composer: &mut StandardComposer,
    state: &[Variable; 3],
) -> [Variable; 3] {
    let montgomery_two = BlsScalar([
        17179869180,
        12756850513266774020,
        3681868479150465002,
        3479420709561305823,
    ]);

    let two = composer.add_witness_to_circuit_description(montgomery_two);

    let x_squared = composer.mul(
        BlsScalar::one(),
        state[0],
        state[0],
        BlsScalar::one(),
        None,
    );

    let var_one = composer.big_mul(
        BlsScalar::one(),
        x_squared,
        x_squared,
        Some((BlsScalar::one(), state[0])),
        BlsScalar::one(),
        None,
    );

    let tuple = composer.big_add(
        (BlsScalar::one(), x_squared),
        (BlsScalar::one(), state[0]),
        Some((BlsScalar::one(), two)),
        BlsScalar::one(),
        None,
    );

    let var_two =
        composer.mul(BlsScalar::one(), state[1], tuple, BlsScalar::one(), None);

    let montgomery_three = BlsScalar([
        25769803770,
        688531696190609414,
        14746174755580473312,
        5219131064341958734,
    ]);

    let three = composer.add_witness_to_circuit_description(montgomery_three);

    let montgomery_four = BlsScalar([
        34359738360,
        7066956952823996424,
        7363736958300930005,
        6958841419122611646,
    ]);

    let four = composer.add_witness_to_circuit_description(montgomery_four);

    // x3 ·(x2 +α2 ·x2 +β2))

    let tuple_one =
        composer.mul(BlsScalar::one(), three, state[1], BlsScalar::one(), None);
    let tuple_two = composer.big_add(
        (BlsScalar::one(), tuple_one),
        (BlsScalar::one(), four),
        Some((BlsScalar::one(), state[1])),
        BlsScalar::one(),
        None,
    );

    let var_three = composer.mul(
        BlsScalar::one(),
        tuple_two,
        state[2],
        BlsScalar::one(),
        None,
    );

    [var_one, var_two, var_three]
}

/// This function computes the in-circuit bars function,
/// as part of the hashing gadget
pub fn bar_gadget(composer: &mut StandardComposer, state: &[Variable; 3]) {
    let zero = composer.add_witness_to_circuit_description(BlsScalar::zero());
}

/// In-circuit concrete function as part of the Zelbet hashing
/// gadget with t = 3 and MDS matrix M = circ(2, 1, 1).
pub fn concrete_gadget(
    composer: &mut StandardComposer,
    state: &[Variable; 3],
    constants: &[Variable; 3],
) -> [Variable; 3] {
    let montgomery_two = BlsScalar([
        17179869180,
        12756850513266774020,
        3681868479150465002,
        3479420709561305823,
    ]);

    let two = composer.add_witness_to_circuit_description(montgomery_two);

    // out0 = 2*u[0] + u[1] + u[2] + c[0];
    let a0 =
        composer.mul(BlsScalar::one(), two, state[0], BlsScalar::one(), None);
    let b0 = composer.big_add(
        (BlsScalar::one(), a0),
        (BlsScalar::one(), state[1]),
        Some((BlsScalar::one(), two)),
        BlsScalar::one(),
        None,
    );
    let c0 = composer.big_add(
        (BlsScalar::one(), b0),
        (BlsScalar::one(), state[2]),
        Some((BlsScalar::one(), two)),
        BlsScalar::one(),
        None,
    );
    let out0 = composer.big_add(
        (BlsScalar::one(), c0),
        (BlsScalar::one(), constants[0]),
        Some((BlsScalar::one(), two)),
        BlsScalar::one(),
        None,
    );

    // out1 = u[0] + 2*u[1] + u[2] + c[1];
    let a1 =
        composer.mul(BlsScalar::one(), two, state[1], BlsScalar::one(), None);
    let b1 = composer.big_add(
        (BlsScalar::one(), a1),
        (BlsScalar::one(), state[0]),
        Some((BlsScalar::one(), two)),
        BlsScalar::one(),
        None,
    );
    let c1 = composer.big_add(
        (BlsScalar::one(), b1),
        (BlsScalar::one(), state[2]),
        Some((BlsScalar::one(), two)),
        BlsScalar::one(),
        None,
    );
    let out1 = composer.big_add(
        (BlsScalar::one(), c1),
        (BlsScalar::one(), constants[1]),
        Some((BlsScalar::one(), two)),
        BlsScalar::one(),
        None,
    );

    // out2 = u[0] + u[1] + 2*u[2] + c[2];
    let a2 =
        composer.mul(BlsScalar::one(), two, state[2], BlsScalar::one(), None);
    let b2 = composer.big_add(
        (BlsScalar::one(), a2),
        (BlsScalar::one(), state[0]),
        Some((BlsScalar::one(), two)),
        BlsScalar::one(),
        None,
    );
    let c2 = composer.big_add(
        (BlsScalar::one(), b2),
        (BlsScalar::one(), state[1]),
        Some((BlsScalar::one(), two)),
        BlsScalar::one(),
        None,
    );
    let out2 = composer.big_add(
        (BlsScalar::one(), c2),
        (BlsScalar::one(), constants[2]),
        Some((BlsScalar::one(), two)),
        BlsScalar::one(),
        None,
    );

    return [out0, out1, out2];
}

#[cfg(test)]
mod tests {
    use super::*;
    use dusk_plonk::constraint_system::StandardComposer;

    #[test]
    fn test_bricks_gadget() {
        let mut composer = StandardComposer::new();
        let one = composer.add_witness_to_circuit_description(BlsScalar::one());
        let output = brick_gagdet(&mut composer, &[one, one, one]);
        composer.constrain_to_constant(
            output[2],
            BlsScalar::from(16384),
            Some(BlsScalar::zero()),
        );
    }
}
