// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use crate::constants::{MONTGOMERY_TWO, S_I_DECOMPOSITION_MONTGOMERY};

pub use super::zelbet::*;
use crate::constants::CONSTANTS_BLS;
use dusk_bytes::*;
use dusk_plonk::{
    constraint_system::{StandardComposer, Variable},
    prelude::BlsScalar,
};

/// Sponge design for in circuit reinforced concrete, with arbitrary length
/// input and output. Input length is read by the function, and output length is
/// an input parameter.
pub fn sponge_zelbet_gadget(
    composer: &mut StandardComposer,
    input: Vec<Variable>,
    length_out: usize,
) -> Vec<Variable> {
    // These constants are needed in reinforced concrete, so they are recorded
    // as variables here for efficiency
    let zero = composer.add_witness_to_circuit_description(BlsScalar::zero());
    let one = composer.add_witness_to_circuit_description(BlsScalar::one());
    let two = composer.add_witness_to_circuit_description(MONTGOMERY_TWO);
    // We don't want to add the constants for the concrete gadget as inputs for
    // each call of zelbet, so we add them once in the sponge
    let mut constants_for_rounds = [one; 18];
    (0..6).for_each(|k| {
        (0..3).for_each(|j| {
            constants_for_rounds[3 * k + j] = composer
                .add_witness_to_circuit_description(CONSTANTS_BLS[k][j]);
        })
    });

    // Initialise the s_i values as variables
    let mut s_i_decomposition = [input[0]; 27];
    (0..27).for_each(|k| {
        s_i_decomposition[k] = composer.add_witness_to_circuit_description(
            S_I_DECOMPOSITION_MONTGOMERY[k],
        );
    });

    // Initialise input as mutable
    let mut input = input.clone();

    // Determine input length to carry out correct number of rounds, and pad if
    // it is not a multiple of two
    let mut length_in = input.len();
    if length_in % 2 == 1 {
        input.push(one);
        length_in = length_in + 1;
    }

    // Fixed starting constant values as defined in reinforced concrete paper
    let mut state = [input[0], input[1], one];

    state = zelbet_gadget(
        composer,
        &state,
        s_i_decomposition,
        constants_for_rounds,
        zero,
        one,
        two,
    );

    // Repeat cycle of adding the two relevant scalars together and then hashing
    (1..length_in / 2).for_each(|k| {
        // Field addition of the two scalars together
        state[0] = composer.big_add(
            (BlsScalar::one(), state[0]),
            (BlsScalar::one(), input[2 * k]),
            None,
            BlsScalar::zero(),
            None,
        );
        state[1] = composer.big_add(
            (BlsScalar::one(), state[1]),
            (BlsScalar::one(), input[2 * k + 1]),
            None,
            BlsScalar::zero(),
            None,
        );
        // Conduct the next round of hashing
        state = zelbet_gadget(
            composer,
            &state,
            s_i_decomposition,
            constants_for_rounds,
            zero,
            one,
            two,
        );
    });

    // Initialise output vector
    let mut output = vec![state[0]; length_out];
    if length_out > 1 {
        output[1] = state[1];
    };

    if length_out % 2 == 0 {
        (1..length_out / 2).for_each(|k| {
            state = zelbet_gadget(
                composer,
                &state,
                s_i_decomposition,
                constants_for_rounds,
                zero,
                one,
                two,
            );
            output[2 * k] = state[0];
            output[2 * k + 1] = state[1];
        })
    } else if length_out % 2 == 1 && length_out > 1 {
        (1..(length_out - 1) / 2).for_each(|k| {
            state = zelbet_gadget(
                composer,
                &state,
                s_i_decomposition,
                constants_for_rounds,
                zero,
                one,
                two,
            );
            output[2 * k] = state[0];
            output[2 * k + 1] = state[1];
        });
        state = zelbet_gadget(
            composer,
            &state,
            s_i_decomposition,
            constants_for_rounds,
            zero,
            one,
            two,
        );
        output[length_out - 1] = state[0];
    }

    output
}

/// Sponge design for out of circuit reinforced concrete, with arbitrary length
/// input and output. Input length is read by the function, and output length is
/// an input parameter.
pub fn sponge_zelbet_out_of(
    input: Vec<BlsScalar>,
    length_out: usize,
) -> Vec<BlsScalar> {
    // Initialise input as mutable
    let mut input = input.clone();

    // Determine input length to carry out correct number of rounds, and pad if
    // it is not a multiple of two
    let mut length_in = input.len();
    if length_in % 2 == 1 {
        input.push(BlsScalar::one());
        length_in = length_in + 1;
    }

    // Fixed starting constant values as defined in reinforced concrete paper
    let mut state = [input[0], input[1], BlsScalar::one()];

    state = zelbet_out_of_circuit(state);

    // Repeat cycle of adding the two relevant scalars together and then hashing
    (1..length_in / 2).for_each(|k| {
        // Add the relevant scalars together
        state[0] = state[0] + input[2 * k];
        state[1] = state[1] + input[2 * k + 1];
        // Conduct the next round of hashing
        state = zelbet_out_of_circuit(state);
    });

    // Initialise output vector
    let mut output = vec![state[0]; length_out];
    if length_out > 1 {
        output[1] = state[1];
    };

    if length_out % 2 == 0 {
        (1..length_out / 2).for_each(|k| {
            state = zelbet_out_of_circuit(state);
            output[2 * k] = state[0];
            output[2 * k + 1] = state[1];
        })
    } else if length_out % 2 == 1 && length_out > 1 {
        (1..(length_out - 1) / 2).for_each(|k| {
            state = zelbet_out_of_circuit(state);
            output[2 * k] = state[0];
            output[2 * k + 1] = state[1];
        });
        state = zelbet_out_of_circuit(state);
        output[length_out - 1] = state[0];
    }

    output
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gadget_tester;
    use dusk_plonk::plookup::PlookupTable4Arity;
    use test::Bencher;

    #[bench]
    fn bench_sponge_out(b: &mut Bencher) {
        let state = vec![BlsScalar::from(1); 5];
        let mut length_out = 5;

        b.iter(|| {
            sponge_zelbet_out_of(state.clone(), length_out);
        });
    }

    #[bench]
    fn bench_sponge_in(b: &mut Bencher) {
        let mut composer = StandardComposer::new();
        let hash_table = PlookupTable4Arity::create_hash_table();
        composer.append_lookup_table(&hash_table);
        let one = composer.add_input(BlsScalar::one());
        let minus_one = composer.add_input(-BlsScalar::one());
        let in3 = composer.add_input(BlsScalar::from(23848872923));
        let in4 = composer.add_input(BlsScalar::from(298375439085));
        let in5 = composer.add_input(-BlsScalar::from(45));
        let input = vec![one, minus_one, in3, in4, in5];

        b.iter(|| {
            sponge_zelbet_gadget(&mut composer, input.clone(), 5);
        });
    }

    // Currently nothing to actually test this result against, this test simply
    // checks whether the function runs or not. Should add results from an
    // independent python programme to compare against
    #[test]
    fn test_sponge_in_circuit() {
        let res = gadget_tester(
            |composer| {
                // This code will be used to rederive hard coded comparison once
                // new conrete constants are added let in_out_of
                // = vec![     BlsScalar::one(),
                //     -BlsScalar::one(),
                //     BlsScalar::from(23848872923),
                //     BlsScalar::from(298375439085),
                //     -BlsScalar::from(45),
                // ];
                // let out_result = sponge_zelbet_out_of(in_out_of, 5);
                // (0..5).for_each(|k| {
                //     println!("BlsScalar({:?}),", out_result[k].0);
                // });
                let hash_table = PlookupTable4Arity::create_hash_table();
                composer.append_lookup_table(&hash_table);
                let one = composer.add_input(BlsScalar::one());
                let minus_one = composer.add_input(-BlsScalar::one());
                let in3 = composer.add_input(BlsScalar::from(23848872923));
                let in4 = composer.add_input(BlsScalar::from(298375439085));
                let in5 = composer.add_input(-BlsScalar::from(45));
                let input = vec![one, minus_one, in3, in4, in5];
                println!("circuit size: {:?}", composer.circuit_size());
                let result = sponge_zelbet_gadget(composer, input, 5);
                println!("circuit size: {:?}", composer.circuit_size());

                // Compare output values to output from out of circuit version
                // to help ensure consistency
                composer.constrain_to_constant(
                    result[0],
                    BlsScalar([
                        3484747639857072351,
                        14140243809318516317,
                        2677010683091861284,
                        5989822391440071232,
                    ]),
                    None,
                );
                composer.constrain_to_constant(
                    result[1],
                    BlsScalar([
                        15789496571931201497,
                        5734252351908561059,
                        14480882572942963540,
                        7231953833215603379,
                    ]),
                    None,
                );
                composer.constrain_to_constant(
                    result[2],
                    BlsScalar([
                        12475486859069551236,
                        14704606153301152199,
                        13339627484282619079,
                        5821826886014145595,
                    ]),
                    None,
                );
                composer.constrain_to_constant(
                    result[3],
                    BlsScalar([
                        15842066148917484615,
                        17957238005947866086,
                        17629428195131741089,
                        861367079814241638,
                    ]),
                    None,
                );
                composer.constrain_to_constant(
                    result[4],
                    BlsScalar([
                        15069830105509539134,
                        11694502824869259129,
                        12399769568792865309,
                        5364747168473928464,
                    ]),
                    None,
                );

                // This code will be used to rederive hard coded comparison once
                // new conrete constants are added let in_out_of
                // = vec![     BlsScalar::from(2),
                //     -BlsScalar::one(),
                //     BlsScalar::from(23848872923),
                //     BlsScalar::from(298375439085),
                // ];
                // let out_result = sponge_zelbet_out_of(in_out_of, 4);
                // (0..5).for_each(|k| {
                //     println!("BlsScalar({:?}),", out_result[k].0);
                // });

                let two = composer.add_input(BlsScalar::from(2));
                println!("circuit size: {:?}", composer.circuit_size());
                let result2 = sponge_zelbet_gadget(
                    composer,
                    vec![two, minus_one, in3, in4],
                    4,
                );
                println!("circuit size: {:?}", composer.circuit_size());

                // Compare output values to output from out of circuit version
                // to help ensure consistency
                composer.constrain_to_constant(
                    result2[0],
                    BlsScalar([
                        50388045448079473,
                        384051309187165497,
                        17074440178406220027,
                        4622546261568619502,
                    ]),
                    None,
                );
                composer.constrain_to_constant(
                    result2[1],
                    BlsScalar([
                        712597510735913230,
                        13604629047046587810,
                        6697632738670923132,
                        2667435796250906486,
                    ]),
                    None,
                );
                composer.constrain_to_constant(
                    result2[2],
                    BlsScalar([
                        4747679144187806666,
                        17057363089776844537,
                        14268833567824355452,
                        7245403559243854770,
                    ]),
                    None,
                );
                composer.constrain_to_constant(
                    result2[3],
                    BlsScalar([
                        15177100986493334512,
                        6414273404242081015,
                        11742282880708228053,
                        2437632187824534994,
                    ]),
                    None,
                );
            },
            5000,
        );
        assert!(res.is_ok());
    }

    // Currently nothing to actually test this result against, this test simply
    // checks whether the function runs or not. Should add results from an
    // independent python programme to compare against
    #[test]
    fn test_sponge_zelbet() {
        let state = vec![BlsScalar::from(1); 5];
        let mut length_out = 5;
        let output = sponge_zelbet_out_of(state, length_out);
        (0..length_out).for_each(|k| {
            println!("value {} is {:?}", k + 1, output[k].0);
        });

        length_out = 4;
        let output2 =
            sponge_zelbet_out_of(vec![BlsScalar::from(2); 4], length_out);
        (0..length_out).for_each(|k| {
            println!("second value {} is {:?}", k + 1, output2[k].0);
        });
    }
}
