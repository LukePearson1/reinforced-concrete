//! Here the constants and fundamental building blocks, which make up the moving parts
//! of the hash function Reinforced Concrete, are defined.

extern crate dusk_bls12_381 as BLS;
use std::convert::TryInto;

use bigint::U256 as u256;
use dusk_bytes::Serializable;
use BLS::BlsScalar as Scalar;

pub mod constants;
use constants::{decomposition_inverses_mont, vu256, SboxBLS};

// Convert representation from tuple in (Z_{s1} x ... x Z_{s_n}) to single element
fn compute_whole_representation(decomposition: [u256; 54]) -> Scalar {
    let mut single: u256 = u256::zero();
    // Note that decomposition[53] is s_1, not s_n, so decomposition[1] is s_n
    for k in (0..27).rev() {
        single = single + decomposition[2 * k];
        if k > 0 {
            single = single * decomposition[2 * k - 1]
        }
    }
    Scalar::from_raw(single.0)
}

// S-box used in bar function
fn small_s_box(x: u256) -> u256 {
    let mut y = u256::zero();
    if x < vu256 {
        let index: u32 = x.try_into().unwrap();
        y = SboxBLS[index as usize];
    } else {
        y = x
    }
    y
}

// Lookup-table-based Sbox
fn bar(state: [Scalar; 3], decomposition: [u256; 54], v: u256) -> [Scalar; 3] {
    let mut new_state = [u256::zero(); 3];
    let mut neww_state = [Scalar::zero(); 3];
    let mut nibbles = [u256::zero(); 27];
    for i in 0..3 {
        // 1. Decomposition
        // Get state value that we are decomposing and manipulating
        let mut intermediate = u256(state[i].reduce().0);
        let mut value = u256::zero();
        for k in 0..27 {
            value = intermediate % decomposition[2 * k + 1];
            // Reduce intermediate representation
            if k < 26 {
                // Convert to BLS scalar form to make use of fast modular multiplication (rather than dividing)
                let intermediate_scalar: Scalar =
                    Scalar::from_raw((intermediate - value).0) * decomposition_inverses_mont[k];
                intermediate = u256(intermediate_scalar.reduce().0);
            } else {
                value = intermediate
            }
            // 2. S-box
            nibbles[k] = small_s_box(value);
        }

        // 3. Composition
        for k in (0..27).rev() {
            new_state[i] = new_state[i] + nibbles[k];
            if k > 0 {
                new_state[i] = new_state[i] * decomposition[2 * k - 1];
            }
        }
        neww_state[i] = Scalar(new_state[i].0);
    }
    neww_state
}

// Element-wise power function
fn brick(state: [Scalar; 3]) -> [Scalar; 3] {
    let mut new_state = [Scalar::zero(); 3];
    let squaring = *Scalar::from(2).internal_repr();
    let x_squared = state[0].pow(&squaring);
    new_state[0] = &x_squared.pow(&squaring) * state[0];
    new_state[1] = state[1] * (&x_squared + state[0] + Scalar::from(2));
    new_state[2] = state[2] * (state[1] * state[1] + Scalar::from(3) * state[1] + Scalar::from(4));
    new_state
}

// Apply affine transformation to state via MDS matrix multiplication
fn concrete(state: [Scalar; 3], matrix: [[u32; 3]; 3], constants: [Scalar; 3]) -> [Scalar; 3] {
    let mut new_state = constants;
    // matrix multiplication
    for i in 0..3 {
        for j in 0..3 {
            new_state[i] += Scalar::from(matrix[i][j] as u64) * state[j];
        }
    }
    new_state
}

// Reinforced concrete hash function, taking in the hash parameters and three-element item
// to be hashed, and outputting the hash value (three BLS scalar elements)
pub fn full_hash(
    state: [Scalar; 3],
    decomposition: [u256; 54],
    v: u256,
    matrix: [[u32; 3]; 3],
    constants: [[Scalar; 3]; 6],
) -> [Scalar; 3] {
    let mut new_state = concrete(state, matrix, constants[0].clone());
    new_state = brick(new_state);
    new_state = concrete(new_state, matrix, constants[1].clone());
    new_state = brick(new_state);
    new_state = concrete(new_state, matrix, constants[2].clone());
    new_state = bar(new_state, decomposition, v);
    new_state = concrete(new_state, matrix, constants[3].clone());
    new_state = brick(new_state);
    new_state = concrete(new_state, matrix, constants[4].clone());
    new_state = brick(new_state);
    new_state = concrete(new_state, matrix, constants[5].clone());
    new_state
}

#[cfg(test)]
mod tests {
    use super::*;
    use constants::{constantsBLS, decompositionBLSS, decomposition_inverses_mont, matrixBLS};

    #[test]
    fn decomposition_inverses_correct() {
        for k in 0..27 {
            let product =
                Scalar::from_raw(decompositionBLSS[2 * k + 1].0) * decomposition_inverses_mont[k];
            assert_eq!(product, Scalar::one());
        }
    }

    #[test]
    fn test1_1() {
        let input = [compute_whole_representation(decompositionBLSS) - Scalar::from(1); 3];
        let output1 = bar(input, decompositionBLSS, vu256);
        assert_eq!(input, output1);
    }

    #[test]
    fn test1_2() {
        let input = [compute_whole_representation(decompositionBLSS) - Scalar::from(1); 3];
        let output2 = brick(input);
        assert_eq!(input, output2);
    }

    #[test]
    fn test1_3() {
        let input = [compute_whole_representation(decompositionBLSS) - Scalar::from(1); 3];
        let output3 = concrete(input, matrixBLS, constantsBLS[1]);
        assert_eq!(input, output3);
    }

    #[test]
    fn test2() {
        let input2 = [compute_whole_representation(decompositionBLSS) - Scalar::from(2); 3];
        let output2 = bar(input2, decompositionBLSS, vu256);
        assert_eq!(input2, output2);
    }

    #[test]
    fn test3() {
        let input3 = [Scalar::from(0); 3];
        let output3 = bar(input3, decompositionBLSS, vu256);
        assert_eq!(input3, output3);
    }

    #[test]
    fn test4() {
        let input4 = [Scalar::from(1); 3];
        let output4 = bar(input4, decompositionBLSS, vu256);
        assert_eq!(input4, output4);
    }

    #[test]
    fn test5() {
        let input5 = [Scalar::from(5); 3];
        let output5 = bar(input5, decompositionBLSS, vu256);
        assert_eq!(input5, output5);
    }
}
