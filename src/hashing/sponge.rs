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
            constants_for_rounds[3 * k + j] =
                composer.add_input(CONSTANTS_BLS[k][j]);
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
            BlsScalar::zero(),
        );
        state[1] = composer.big_add(
            (BlsScalar::one(), state[1]),
            (BlsScalar::one(), input[2 * k + 1]),
            None,
            BlsScalar::zero(),
            BlsScalar::zero(),
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

/// Sponge design for reinforced concrete, with arbitrary length input and
/// output. Input length is read by the function, and output length is an input
/// parameter.
pub fn sponge_zelbet(
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

    // Currently nothing to actually test this result against, this test simply
    // checks whether the function runs or not. Should add results from an
    // independent python programme to compare against
    #[test]
    fn test_sponge_gadget() {
        let res = gadget_tester(
            |composer| {
                let hash_table = PlookupTable4Arity::create_hash_table();
                composer.append_lookup_table(&hash_table);
                let one = composer.add_input(BlsScalar::one());
                let result = sponge_zelbet_gadget(composer, vec![one; 5], 5);

                // Compare output values to output from out of circuit version
                // to help ensure consistency
                composer.constrain_to_constant(
                    result[0],
                    BlsScalar([
                        2948929853694419607,
                        7328611936786299993,
                        9919458621989105704,
                        6130852493708878965,
                    ]),
                    BlsScalar::zero(),
                );
                composer.constrain_to_constant(
                    result[1],
                    BlsScalar([
                        17299957155812118625,
                        17339778072127637790,
                        11924609608675978405,
                        2154428716724127434,
                    ]),
                    BlsScalar::zero(),
                );
                composer.constrain_to_constant(
                    result[2],
                    BlsScalar([
                        8071993871493204166,
                        11515087050602589091,
                        18009909627082794530,
                        7421486182814624372,
                    ]),
                    BlsScalar::zero(),
                );
                composer.constrain_to_constant(
                    result[3],
                    BlsScalar([
                        16302036423968262150,
                        41607565966230019,
                        17585612505115477194,
                        7596324467329926765,
                    ]),
                    BlsScalar::zero(),
                );
                composer.constrain_to_constant(
                    result[4],
                    BlsScalar([
                        3279543462538777620,
                        16893969979408605182,
                        14803853461423869094,
                        1742924770726305161,
                    ]),
                    BlsScalar::zero(),
                );

                let two = composer.add_input(BlsScalar::from(2));
                let result2 = sponge_zelbet_gadget(composer, vec![two; 4], 4);

                // Compare output values to output from out of circuit version
                // to help ensure consistency
                composer.constrain_to_constant(
                    result2[0],
                    BlsScalar([
                        11579792014922627373,
                        2934221241967416843,
                        8345478288467627804,
                        7623670637766636891,
                    ]),
                    BlsScalar::zero(),
                );
                composer.constrain_to_constant(
                    result2[1],
                    BlsScalar([
                        4874119729451523912,
                        8725796636719379660,
                        7985091090312449276,
                        8228914025318238834,
                    ]),
                    BlsScalar::zero(),
                );
                composer.constrain_to_constant(
                    result2[2],
                    BlsScalar([
                        14622228267955924779,
                        9615659355454524416,
                        2250804130334070584,
                        6587213602798012345,
                    ]),
                    BlsScalar::zero(),
                );
                composer.constrain_to_constant(
                    result2[3],
                    BlsScalar([
                        18368122655394568764,
                        410839651173137500,
                        18398795598074206863,
                        2720386774260300492,
                    ]),
                    BlsScalar::zero(),
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
        let output = sponge_zelbet(state, length_out);
        (0..length_out).for_each(|k| {
            println!("value {} is {:?}", k + 1, output[k].0);
        });

        length_out = 4;
        let output2 = sponge_zelbet(vec![BlsScalar::from(2); 4], length_out);
        (0..length_out).for_each(|k| {
            println!("second value {} is {:?}", k + 1, output2[k].0);
        });
    }
}
