// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! This file contains the circuit implementation of the
//! zelbet hash function

<<<<<<< HEAD
use crate::constants::DECOMPOSITION_S_I;
=======
use crate::constants::{DECOMPOSITION_S_I, INVERSES_S_I, MONTGOMERY_TWO};
use bigint::U256 as u256;
>>>>>>> acc8b02e8ff9a706e20e0d23ff8e1a2461ec998d
use dusk_plonk::constraint_system::{StandardComposer, Variable};
use dusk_plonk::prelude::*;

/// This function computes the in-circuit brick function,
/// as part of the hashing gadget
pub fn brick_gadget(
    composer: &mut StandardComposer,
    state: &[Variable; 3],
    two: Variable,
) -> [Variable; 3] {
    let x_squared = composer.mul(
        BlsScalar::one(),
        state[0],
        state[0],
        BlsScalar::zero(),
        BlsScalar::zero(),
    );
    let x_fourth = composer.mul(
        BlsScalar::one(),
        x_squared,
        x_squared,
        BlsScalar::zero(),
        BlsScalar::zero(),
    );
    let var_one = composer.mul(
        BlsScalar::one(),
        x_fourth,
        state[0],
        BlsScalar::zero(),
        BlsScalar::zero(),
    );

    let tuple = composer.big_add(
        (BlsScalar::one(), x_squared),
        (BlsScalar::one(), state[0]),
        Some((BlsScalar::one(), two)),
        BlsScalar::zero(),
        BlsScalar::zero(),
    );
    let var_two = composer.mul(
        BlsScalar::one(),
        state[1],
        tuple,
        BlsScalar::zero(),
        BlsScalar::zero(),
    );

    // x3 ·(x2^2 +α2 ·x2 +β2))
    let y_squared_plus_4 = composer.mul(
        BlsScalar::one(),
        state[1],
        state[1],
        BlsScalar::from(4),
        BlsScalar::zero(),
    );
    let tuple_one = composer.big_add(
        (BlsScalar::one(), y_squared_plus_4),
        (BlsScalar::from(3), state[1]),
        None,
        BlsScalar::zero(),
        BlsScalar::zero(),
    );
    let var_three = composer.mul(
        BlsScalar::one(),
        tuple_one,
        state[2],
        BlsScalar::zero(),
        BlsScalar::zero(),
    );

    [var_one, var_two, var_three]
}

/// In-circuit concrete function as part of the Zelbet hashing
/// gadget with t = 3 and MDS matrix M = circ(2, 1, 1).
pub fn concrete_gadget(
    composer: &mut StandardComposer,
    state: &[Variable; 3],
    constants: &[Variable; 3],
) -> [Variable; 3] {
    // out0 = 2*u[0] + u[1] + u[2] + c[0];
    let a0 = composer.big_add(
        (BlsScalar::from(2), state[0]),
        (BlsScalar::one(), state[1]),
        Some((BlsScalar::one(), state[2])),
        BlsScalar::zero(),
        BlsScalar::zero(),
    );
    let out0 = composer.add(
        (BlsScalar::one(), a0),
        (BlsScalar::one(), constants[0]),
        BlsScalar::zero(),
        BlsScalar::zero(),
    );

    // out1 = u[0] + 2*u[1] + u[2] + c[1];
    let a1 = composer.big_add(
        (BlsScalar::one(), state[0]),
        (BlsScalar::from(2), state[1]),
        Some((BlsScalar::one(), state[2])),
        BlsScalar::zero(),
        BlsScalar::zero(),
    );
    let out1 = composer.add(
        (BlsScalar::one(), a1),
        (BlsScalar::one(), constants[1]),
        BlsScalar::zero(),
        BlsScalar::zero(),
    );

    // out2 = u[0] + u[1] + 2*u[2] + c[2];
    let a2 = composer.big_add(
        (BlsScalar::one(), state[0]),
        (BlsScalar::one(), state[1]),
        Some((BlsScalar::from(2), state[2])),
        BlsScalar::zero(),
        BlsScalar::zero(),
    );
    let out2 = composer.add(
        (BlsScalar::one(), a2),
        (BlsScalar::one(), constants[2]),
        BlsScalar::zero(),
        BlsScalar::zero(),
    );

    [out0, out1, out2]
}

<<<<<<< HEAD
/// Bar function
pub fn bar_gadget(
    composer: &mut StandardComposer,
    input: Variable,
    s_i_decomposition: [Variable; 27],
    one: Variable,
    two: Variable,
) -> Variable {
    let mut tuple = composer.decomposition_gadget(input, s_i_decomposition);

    // Initialise the constraints
    let mut c_i = [input; 27];
    let mut conditional = false;
    let mut z_i = [input; 27];
    (0..27).for_each(|k| {
        let result = composer.s_box_and_constraints(
            tuple[k],
            (27 - k) as u64,
            conditional,
            one,
            two,
        );
        tuple[k] = result.0;
        c_i[k] = result.1;
        conditional = result.2;
        z_i[k] = result.3;
    });

    // Constraint checks for c_i, bearing in mind that c_i[0] = c_27
    (5..=26).step_by(3).into_iter().for_each(|k| {
        composer.plookup_gate(
            c_i[k],
            c_i[k - 1],
            c_i[k - 2],
            Some(c_i[k - 3]),
            BlsScalar::zero(),
        );
    });
    composer.plookup_gate(
        c_i[3],
        c_i[2],
        c_i[1],
        Some(c_i[0]),
        BlsScalar::zero(),
    );

    // Constraint checks for z_i, bearing in mind that z_i[0] = z_27
    (6..=26).step_by(4).for_each(|k| {
        composer.plookup_gate(
            z_i[k],
            z_i[k - 1],
            z_i[k - 2],
            Some(z_i[k - 3]),
            BlsScalar::zero(),
        );
    });
    composer.plookup_gate(
        z_i[3],
        z_i[2],
        z_i[1],
        Some(z_i[0]),
        BlsScalar::zero(),
    );

    let mut accumulator_var = composer.add_input(BlsScalar::zero());
    accumulator_var = composer.big_add(
        (BlsScalar::one(), accumulator_var),
        (BlsScalar::one(), tuple[26]),
        None,
        BlsScalar::zero(),
        BlsScalar::zero(),
    );
    (1..=26).rev().for_each(|k| {
        let s_i_var =
            composer.add_input(BlsScalar::from_raw(DECOMPOSITION_S_I[k - 1].0));
        accumulator_var = composer.big_mul(
            BlsScalar::one(),
            accumulator_var,
            s_i_var,
            Some((BlsScalar::one(), tuple[k - 1])),
            BlsScalar::zero(),
            BlsScalar::zero(),
        );
    });

    accumulator_var
}
=======
// /// Bar function
// pub fn bar_gadget(composer: &mut StandardComposer, input: Variable) -> Variable {
//     let mut tuple = composer.decomposition_gadget(input, DECOMPOSITION_S_I, INVERSES_S_I);

//     // let s_box_table = PlookupTable4Arity::s_box_table();
//     // composer.lookup_table = s_box_table;
//     (0..27).for_each(|k| {
//         tuple[k] = composer.s_box(tuple[k]);
//     });

//     let mut accumulator_var = composer.add_input(BlsScalar::zero());
//     (1..27).rev().for_each(|k| {
//         if k == 26 {
//             accumulator_var = composer.big_add(
//                 (BlsScalar::one(), accumulator_var),
//                 (BlsScalar::one(), tuple[k]),
//                 None,
//                 BlsScalar::zero(),
//                 None,
//             );
//         }
//         let s_i_var = composer.add_input(BlsScalar::from_raw(DECOMPOSITION_S_I[k - 1].0));
//         accumulator_var = composer.big_mul(
//             BlsScalar::one(),
//             accumulator_var,
//             s_i_var,
//             Some((BlsScalar::one(), tuple[k - 1])),
//             BlsScalar::zero(),
//             None,
//         );
//     });

//     accumulator_var
// }
>>>>>>> acc8b02e8ff9a706e20e0d23ff8e1a2461ec998d

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{constants::MATRIX_BLS, hashing::zelbet::brick};
    use crate::{gadget_tester, hashing::zelbet::concrete};
    use dusk_plonk::plookup::PlookupTable4Arity;

    #[test]
    fn test_bar() {
        let res = gadget_tester(
            |composer| {
                let hash_table = PlookupTable4Arity::create_hash_table();
                composer.append_lookup_table(&hash_table);
                let one = composer.add_input(BlsScalar::one());
                let two = composer.add_input(BlsScalar::from(2));
                let mut s_i_decomposition = [one; 27];
                (0..27).for_each(|k| {
                    s_i_decomposition[k] =
                        composer.add_input(DECOMPOSITION_S_I[k]);
                });
                // Check that the output is what we expected (in Montgomery)
                let output =
                    bar_gadget(composer, one, s_i_decomposition, one, two);
                let expected = BlsScalar([
                    2921300856332839541,
                    8943181998193365483,
                    12554333934768435622,
                    1625679107374292725,
                ]);
                composer.constrain_to_constant(
                    output,
                    expected,
                    BlsScalar::zero(),
                );
                let one_eight_seven = composer.add_input(BlsScalar::from(187));
                let zero = composer.add_input(BlsScalar::from(0));

                // Plookup is designed to not pass if the number of plookup
                // checks is much smaller than the size of the lookup table, so
                // these extra checks are to increase the number of plookup
                // gates
                (0..550).for_each(|_| {
                    composer.plookup_gate(
                        one,
                        zero,
                        one_eight_seven,
                        Some(one),
                        BlsScalar::zero(),
                    );
                    composer.plookup_gate(
                        one,
                        one,
                        one,
                        Some(one),
                        BlsScalar::zero(),
                    );
                });
            },
            3000,
        );
        assert!(res.is_ok());
    }

    #[test]
    fn test_bricks_gadget() {
        let res = gadget_tester(
            |composer| {
                let two = composer
                    .add_witness_to_circuit_description(BlsScalar::from(2));
                let three = composer
                    .add_witness_to_circuit_description(BlsScalar::from(3));
                let four = composer
                    .add_witness_to_circuit_description(BlsScalar::from(4));
                let output = brick_gadget(composer, &[two, three, four], two);

                let output_1 = brick([
                    BlsScalar::from(2),
                    BlsScalar::from(3),
                    BlsScalar::from(4),
                ]);
                // Check in circuit result against out of circuit result
                for i in 0..3 {
                    composer.constrain_to_constant(
                        output[i],
                        output_1[i],
                        BlsScalar::zero(),
                    );
                }

                // Check in circuit result against expected result
                let expected_result = [
                    BlsScalar::from(32),
                    BlsScalar::from(24),
                    BlsScalar::from(88),
                ];
                for i in 0..3 {
                    composer.constrain_to_constant(
                        output[i],
                        expected_result[i],
                        BlsScalar::zero(),
                    );
                }
            },
            32,
        );
        assert!(res.is_ok());
    }

    #[test]
    fn test_concrete_gadget() {
        let res = gadget_tester(
            |composer| {
                let one = composer
                    .add_witness_to_circuit_description(BlsScalar::one());
                let two = composer
                    .add_witness_to_circuit_description(BlsScalar::from(2));
                let three = composer
                    .add_witness_to_circuit_description(BlsScalar::from(3));
                let output =
                    concrete_gadget(composer, &[one, two, three], &[two; 3]);
                let output_1 = concrete(
                    [BlsScalar::one(), BlsScalar::from(2), BlsScalar::from(3)],
                    MATRIX_BLS,
                    [BlsScalar::from(2); 3],
                );

                // Check in circuit result against out of circuit result
                for i in 0..3 {
                    composer.constrain_to_constant(
                        output[i],
                        output_1[i],
                        BlsScalar::zero(),
                    );
                }

                // Check in circuit result against expected result
                let expected_result = [
                    BlsScalar::from(9),
                    BlsScalar::from(10),
                    BlsScalar::from(11),
                ];
                for i in 0..3 {
                    composer.constrain_to_constant(
                        output[i],
                        expected_result[i],
                        BlsScalar::zero(),
                    );
                }
            },
            32,
        );
        assert!(res.is_ok());
    }

    // #[test]
    // fn test_bar_gadget() {
    //     let res = gadget_tester(
    //         |composer| {
    //             let one = composer.add_input(BlsScalar::one());
    //             // Check that the output is what we expected (in Montgomery)
    //             let output = bar_gadget(composer, one);
    //             let expected = BlsScalar([
    //                 2921300856332839541,
    //                 8943181998193365483,
    //                 12554333934768435622,
    //                 1625679107374292725,
    //             ]);
    //             composer.constrain_to_constant(output, expected, BlsScalar::zero());
    //             println!("circuit is {:?}", composer.circuit_size());
    //         },
    //         800,
    //     );
    //     assert!(res.is_ok());
    // }

}
