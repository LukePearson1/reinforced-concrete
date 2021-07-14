// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use dusk_plonk::prelude::*;
use dusk_plonk::proof_system::{Prover, Verifier};
// use rand_core::{Error, OsRng};
use anyhow::{Error, Result};

// Takes a generic gadget function with no auxillary input and
// tests whether it passes an end-to-end test
// pub(crate) fn gadget_tester(
//     gadget: fn(composer: &mut StandardComposer),
//     n: usize,
// ) -> Result<(), Error> {
//     // Common View
//     let public_parameters = PublicParameters::setup(2 * n, &mut OsRng)?;
//     // Provers View
//     let (proof, public_inputs) = {
//         // Create a prover struct
//         let mut prover = Prover::new(b"demo");

//         // Additionally key the transcript
//         prover.key_transcript(b"key", b"additional seed information");

//         // Add gadgets
//         gadget(&mut prover.mut_cs());

//         // Commit Key
//         let (ck, _) = public_parameters
//             .trim(2 * prover.mut_cs().circuit_size().next_power_of_two())?;

//         // Preprocess circuit
//         prover.preprocess(&ck)?;

//         // Once the prove method is called, the public inputs are cleared
//         // So pre-fetch these before calling Prove
//         let public_inputs = prover.mut_cs().construct_dense_pi_vec();

//         // Compute Proof
//         (prover.prove(&ck)?, public_inputs)
//     };
//     // Verifiers view
//     //
//     // Create a Verifier object
//     let mut verifier = Verifier::new(b"demo");

//     // Additionally key the transcript
//     verifier.key_transcript(b"key", b"additional seed information");

//     // Add gadgets
//     gadget(&mut verifier.mut_cs());

//     // Compute Commit and Verifier Key
//     let (ck, vk) = public_parameters
//         .trim(verifier.mut_cs().circuit_size().next_power_of_two())?;

//     // Preprocess circuit
//     verifier.preprocess(&ck)?;

//     // Verify proof
//     verifier.verify(&proof, &vk, &public_inputs)
// }

/// Takes a generic gadget function with no auxillary input and
/// tests whether it passes an end-to-end test
pub(crate) fn gadget_tester(
    gadget: fn(composer: &mut StandardComposer),
    n: usize,
) -> Result<(), Error> {
    // Common View
    let public_parameters =
        PublicParameters::setup(2 * n, &mut rand::thread_rng())?;
    // Provers View
    let (proof, public_inputs, lookup_table) = {
        // Create a prover struct
        let mut prover = Prover::new(b"demo");

        // Additionally key the transcript
        prover.key_transcript(b"key", b"additional seed information");

        // Add gadgets
        gadget(&mut prover.mut_cs());

        // Commit Key
        let (ck, _) = public_parameters
            .trim(2 * prover.mut_cs().circuit_size().next_power_of_two())?;

        // Preprocess circuit
        prover.preprocess(&ck)?;

        // Once the prove method is called, the public inputs are cleared
        // So pre-fetch these before calling Prove
        let public_inputs = prover.mut_cs().public_inputs.clone();
        let lookup_table = prover.mut_cs().lookup_table.clone();

        // Compute Proof
        (prover.prove(&ck)?, public_inputs, lookup_table)
    };
    // Verifiers view
    //
    // Create a Verifier object
    let mut verifier = Verifier::new(b"demo");

    // Additionally key the transcript
    verifier.key_transcript(b"key", b"additional seed information");

    // Add gadgets
    gadget(&mut verifier.mut_cs());

    // Compute Commit and Verifier Key
    let (ck, vk) = public_parameters
        .trim(verifier.mut_cs().circuit_size().next_power_of_two())?;

    // Preprocess circuit
    verifier.preprocess(&ck)?;

    // Verify proof
    verifier.verify(&proof, &vk, &public_inputs, &lookup_table)
}
