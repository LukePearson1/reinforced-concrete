// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! This file contains the circuit implementation of the
//! zelbet hash function
use super::divide_w_recip;
use crate::constants::{DECOMPOSITION_S_I, SBOX_MONTGOMERY};
use bigint::U256 as u256;
use dusk_plonk::constraint_system::{StandardComposer, Variable};
use dusk_plonk::plookup::table::hash_tables::BLS_SCALAR_REAL;
use dusk_plonk::prelude::*;

/// This function computes the in-circuit brick function,
/// as part of the hashing gadget
pub fn brick_gadget(
    composer: &mut StandardComposer,
    state: &[Variable; 3],
    two: Variable,
) -> [Variable; 3] {
    // Finding y_1
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
    let y_1 = composer.mul(
        BlsScalar::one(),
        x_fourth,
        state[0],
        BlsScalar::zero(),
        BlsScalar::zero(),
    );

    // Finding y_2
    let tuple = composer.big_add(
        (BlsScalar::one(), x_squared),
        (BlsScalar::one(), state[0]),
        Some((BlsScalar::one(), two)),
        BlsScalar::zero(),
        BlsScalar::zero(),
    );
    let y_2 = composer.mul(
        BlsScalar::one(),
        state[1],
        tuple,
        BlsScalar::zero(),
        BlsScalar::zero(),
    );

    // Finding y_3
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
    let y_3 = composer.mul(
        BlsScalar::one(),
        tuple_one,
        state[2],
        BlsScalar::zero(),
        BlsScalar::zero(),
    );

    [y_1, y_2, y_3]
}

/// In-circuit concrete function as part of the Zelbet hashing
/// gadget with t = 3 and MDS matrix M = circ(2, 1, 1).
pub fn concrete_gadget(
    composer: &mut StandardComposer,
    state: &[Variable; 3],
    constants: &[Variable; 3],
) -> [Variable; 3] {
    // y_1 = 2*u[0] + u[1] + u[2] + c[0];
    let a0 = composer.big_add(
        (BlsScalar::from(2), state[0]),
        (BlsScalar::one(), state[1]),
        Some((BlsScalar::one(), state[2])),
        BlsScalar::zero(),
        BlsScalar::zero(),
    );
    let y_1 = composer.add(
        (BlsScalar::one(), a0),
        (BlsScalar::one(), constants[0]),
        BlsScalar::zero(),
        BlsScalar::zero(),
    );

    // y_2 = u[0] + 2*u[1] + u[2] + c[1];
    let a1 = composer.big_add(
        (BlsScalar::one(), state[0]),
        (BlsScalar::from(2), state[1]),
        Some((BlsScalar::one(), state[2])),
        BlsScalar::zero(),
        BlsScalar::zero(),
    );
    let y_2 = composer.add(
        (BlsScalar::one(), a1),
        (BlsScalar::one(), constants[1]),
        BlsScalar::zero(),
        BlsScalar::zero(),
    );

    // y_3 = u[0] + u[1] + 2*u[2] + c[2];
    let a2 = composer.big_add(
        (BlsScalar::one(), state[0]),
        (BlsScalar::one(), state[1]),
        Some((BlsScalar::from(2), state[2])),
        BlsScalar::zero(),
        BlsScalar::zero(),
    );
    let y_3 = composer.add(
        (BlsScalar::one(), a2),
        (BlsScalar::one(), constants[2]),
        BlsScalar::zero(),
        BlsScalar::zero(),
    );

    [y_1, y_2, y_3]
}

/// In circuit bar function, making use of decomposition gadget
/// that is defined in PLONK repo
pub fn bar_gadget(
    composer: &mut StandardComposer,
    input: Variable,
    s_i_decomposition: [Variable; 27],
    zero: Variable,
    one: Variable,
    two: Variable,
) -> Variable {
    // Decomposition
    let (mut tuple_mont, tuple_reduced) =
        composer.decomposition_gadget(input, s_i_decomposition);

    // Initialise the constraints
    let mut c_i = [input; 27];
    let mut conditional = false;
    let mut z_i = [input; 27];
    // Calculate s-box permutation, and constraint values. It is important that
    // tuple[26] is done first as this is x_1, and c_1 should be calculated
    // first, not c_27
    (0..27).rev().for_each(|k| {
        let result = s_box_and_constraints(
            composer,
            tuple_mont[k],
            tuple_reduced[k],
            (27 - k) as u64,
            conditional,
            zero,
            one,
            two,
        );
        tuple_mont[k] = result.0;
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

    // Initialise accumulator - this is needed to check constraint number 17
    // from the reinforced concrete paper, that the decomposition's composition
    // does equal the input, x
    let mut accumulator_var = composer.add_input(BlsScalar::zero());
    accumulator_var = composer.big_add(
        (BlsScalar::one(), accumulator_var),
        (BlsScalar::one(), tuple_mont[26]),
        None,
        BlsScalar::zero(),
        BlsScalar::zero(),
    );
    (1..=26).rev().for_each(|k| {
        accumulator_var = composer.big_mul(
            BlsScalar::one(),
            accumulator_var,
            s_i_decomposition[k - 1],
            Some((BlsScalar::one(), tuple_mont[k - 1])),
            BlsScalar::zero(),
            BlsScalar::zero(),
        );
    });

    accumulator_var
}

/// S-box using hash tables, and outputs constraints c_i, z_i and a boolean
/// counter to help determine the c_i. (y_i, c_i, conditional, z_i)
pub fn s_box_and_constraints(
    composer: &mut StandardComposer,
    input_mont: Variable,
    input_reduced: u256,
    counter: u64,
    conditional: bool,
    zero: Variable,
    one: Variable,
    two: Variable,
) -> (Variable, Variable, bool, Variable) {
    // Need to convert input scalar value to non-Montgomery
    // to allow size comparison
    let mut y_i = input_mont;
    let mut c_i = one;
    let mut conditional_new = conditional;
    let mut z_i = zero;
    let mut z_i_val: u64 = 0;
    if input_reduced.0[0] < 659 {
        y_i = composer.add_input(SBOX_MONTGOMERY[input_reduced.0[0] as usize]);
        conditional_new = true;
    } else {
        y_i = input_mont;
        z_i = one;
        z_i_val = 1;
        if input_reduced.0[0] > BLS_SCALAR_REAL[27 - counter as usize].0[0] {
            c_i = two;
            conditional_new = true
        } else if input_reduced.0[0]
            == BLS_SCALAR_REAL[27 - counter as usize].0[0]
        {
            if conditional == true {
                c_i = two;
                conditional_new = true
            } else {
                c_i = zero
            }
        }
    }

    let scaled_z_i = composer.add_input(BlsScalar::from(counter * z_i_val));
    composer.plookup_gate(
        input_mont,
        scaled_z_i,
        y_i,
        Some(c_i),
        BlsScalar::zero(),
    );

    (y_i, c_i, conditional_new, z_i)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::S_I_DECOMPOSITION_MONTGOMERY;
    use crate::hashing::zelbet::brick;
    use crate::{gadget_tester, hashing::zelbet::concrete};
    use dusk_plonk::plookup::PlookupTable4Arity;
    use test::Bencher;

    #[bench]
    fn bench_bar_gadget(b: &mut Bencher) {
        let mut composer = StandardComposer::new();
        let hash_table = PlookupTable4Arity::create_hash_table();
        composer.append_lookup_table(&hash_table);
        let zero = composer.add_input(BlsScalar::zero());
        let one = composer.add_input(BlsScalar::one());
        let two = composer.add_input(BlsScalar::from(2));
        let minus_one = composer.add_input(-BlsScalar::one());
        let mut s_i_decomposition = [one; 27];
        (0..27).for_each(|k| {
            s_i_decomposition[k] = composer.add_witness_to_circuit_description(
                S_I_DECOMPOSITION_MONTGOMERY[k],
            );
        });
        b.iter(|| {
            (0..3).for_each(|_| {
                bar_gadget(
                    &mut composer,
                    minus_one,
                    s_i_decomposition,
                    zero,
                    one,
                    two,
                );
            })
        });
    }

    #[bench]
    fn bench_decomp(b: &mut Bencher) {
        let mut composer = StandardComposer::new();
        let one = composer.add_input(BlsScalar::one());
        let minus_one = composer.add_input(-BlsScalar::one());
        let mut s_i_decomposition = [one; 27];
        (0..27).for_each(|k| {
            s_i_decomposition[k] = composer.add_witness_to_circuit_description(
                S_I_DECOMPOSITION_MONTGOMERY[k],
            );
        });
        b.iter(|| {
            (0..3).for_each(|_| {
                (0..3).for_each(|_| {
                    composer.decomposition_gadget(minus_one, s_i_decomposition);
                });
            })
        });
    }

    #[test]
    fn test_bar_gadget() {
        let res = gadget_tester(
            |composer| {
                let hash_table = PlookupTable4Arity::create_hash_table();
                composer.append_lookup_table(&hash_table);
                let zero = composer.add_input(BlsScalar::zero());
                let one = composer.add_input(BlsScalar::one());
                let two = composer.add_input(BlsScalar::from(2));
                let mut s_i_decomposition = [one; 27];
                (0..27).for_each(|k| {
                    s_i_decomposition[k] =
                        composer.add_input(S_I_DECOMPOSITION_MONTGOMERY[k]);
                });
                // Check bar funciton on input of 1
                let output = bar_gadget(
                    composer,
                    one,
                    s_i_decomposition,
                    zero,
                    one,
                    two,
                );
                let expected = BlsScalar([
                    2921300856332839541,
                    8943181998193365483,
                    12554333934768435622,
                    1625679107374292725,
                ]);
                // Check that the output is what we expected (in Montgomery)
                composer.constrain_to_constant(
                    output,
                    expected,
                    BlsScalar::zero(),
                );

                // Check bar function on input of -5
                let minus_five = composer.add_input(-BlsScalar::from(5));
                let output2 = bar_gadget(
                    composer,
                    minus_five,
                    s_i_decomposition,
                    zero,
                    one,
                    two,
                );
                composer.constrain_to_constant(
                    output2,
                    BlsScalar([
                        18446742991377793276,
                        7975156907413507843,
                        7849958771744875838,
                        2157424182152352530,
                    ]),
                    BlsScalar::zero(),
                );

                // Plookup is designed to not pass if the number of plookup
                // checks is much smaller than the size of the lookup table, so
                // these extra checks are to increase the number of plookup
                // gates
                let one_eight_seven = composer.add_input(BlsScalar::from(187));
                let zero = composer.add_input(BlsScalar::from(0));
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

    #[test]
    fn test_s_box_and_constraints() {
        let res = gadget_tester(
            |composer| {
                let hash_table = PlookupTable4Arity::create_hash_table();
                composer.append_lookup_table(&hash_table);
                let seven_hundred = composer.add_input(BlsScalar::from(700));
                let zero = composer.add_input(BlsScalar::zero());
                let one = composer.add_input(BlsScalar::one());
                let two = composer.add_input(BlsScalar::from(2));
                let prime = composer.add_input(BlsScalar::from(659));
                let counter: u64 = 1;
                let counter2: u64 = 2;
                let conditional = true;
                let output_700 = s_box_and_constraints(
                    composer,
                    seven_hundred,
                    u256::from(700),
                    counter2,
                    conditional,
                    zero,
                    one,
                    two,
                );
                let output_one = s_box_and_constraints(
                    composer,
                    one,
                    u256::from(1),
                    counter,
                    conditional,
                    zero,
                    one,
                    two,
                );
                let output_prime = s_box_and_constraints(
                    composer,
                    prime,
                    u256::from(659),
                    counter2,
                    conditional,
                    zero,
                    one,
                    two,
                );
                let output_prime_false = s_box_and_constraints(
                    composer,
                    prime,
                    u256::from(659),
                    counter,
                    false,
                    zero,
                    one,
                    two,
                );

                // Check that the s-box works as expected
                composer.constrain_to_constant(
                    output_700.0,
                    BlsScalar::from_raw([700, 0, 0, 0]),
                    BlsScalar::zero(),
                );
                composer.constrain_to_constant(
                    output_one.0,
                    BlsScalar::from_raw([187, 0, 0, 0]),
                    BlsScalar::zero(),
                );
                composer.constrain_to_constant(
                    output_prime.0,
                    BlsScalar::from_raw([659, 0, 0, 0]),
                    BlsScalar::zero(),
                );

                (0..1100).for_each(|k| {
                    composer.plookup_gate(
                        prime,
                        one,
                        prime,
                        Some(one),
                        BlsScalar::zero(),
                    );
                });

                // Check that the c_i are output as expected
                composer.constrain_to_constant(
                    output_700.1,
                    BlsScalar::from(2),
                    BlsScalar::zero(),
                );
                composer.constrain_to_constant(
                    output_one.1,
                    BlsScalar::one(),
                    BlsScalar::zero(),
                );
                composer.constrain_to_constant(
                    output_prime.1,
                    BlsScalar::one(),
                    BlsScalar::zero(),
                );
                composer.constrain_to_constant(
                    output_prime_false.1,
                    BlsScalar::one(),
                    BlsScalar::zero(),
                );

                // Check that the counter is output correctly
                assert!(output_700.2);
                assert!(output_one.2);
                assert!(output_prime.2);
                assert!(!output_prime_false.2);

                // Check that z_i is output correctly
                composer.constrain_to_constant(
                    output_700.3,
                    BlsScalar::one(),
                    BlsScalar::zero(),
                );
                composer.constrain_to_constant(
                    output_one.3,
                    BlsScalar::zero(),
                    BlsScalar::zero(),
                );
                composer.constrain_to_constant(
                    output_prime.3,
                    BlsScalar::one(),
                    BlsScalar::zero(),
                );
                composer.constrain_to_constant(
                    output_prime_false.3,
                    BlsScalar::one(),
                    BlsScalar::zero(),
                );
            },
            4000,
        );
        assert!(res.is_ok());
    }

    #[test]
    fn test_s_box_and_constraints_fails() {
        let res = gadget_tester(
            |composer| {
                let hash_table = PlookupTable4Arity::create_hash_table();
                composer.append_lookup_table(&hash_table);
                let one_hundred = composer.add_input(BlsScalar::from(100));
                let zero = composer.add_input(BlsScalar::zero());
                let one = composer.add_input(BlsScalar::one());
                let two = composer.add_input(BlsScalar::from(2));
                let counter: u64 = 1;
                let conditional = true;
                let output = s_box_and_constraints(
                    composer,
                    one_hundred,
                    u256::from(100),
                    counter,
                    conditional,
                    zero,
                    one,
                    two,
                );
                composer.constrain_to_constant(
                    output.0,
                    BlsScalar::from_raw([200, 0, 0, 0]),
                    BlsScalar::zero(),
                );
                composer.constrain_to_constant(
                    output.0,
                    BlsScalar::from_raw([200, 0, 0, 0]),
                    BlsScalar::zero(),
                );
                composer.constrain_to_constant(
                    output.0,
                    BlsScalar::from_raw([200, 0, 0, 0]),
                    BlsScalar::zero(),
                );

                let prime = composer.add_input(BlsScalar::from(659));
                (0..1100).for_each(|k| {
                    composer.plookup_gate(
                        prime,
                        one_hundred,
                        prime,
                        Some(one_hundred),
                        BlsScalar::zero(),
                    );
                });
            },
            2000,
        );
        assert!(res.is_err());
    }
}
