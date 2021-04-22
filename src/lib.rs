//! Here the constants and fundamental building blocks, which make up the moving parts
//! of the hash function Reinforced Concrete, are defined.

// use num_bigint::{BigInt, Sign};
extern crate num_bigint;
use num_bigint::BigUint;

// decomposition = [v_n, s_n, v_{n-1}, s_{n-1} ..., v_1, s_1]
const decompositionBLS: [u32; 54] = [
    15, 688, 17, 699, 0, 692, 1, 709, 12, 708, 4, 677, 13, 685, 21, 693, 5, 673, 11, 703, 28, 701,
    16, 701, 20, 708, 22, 690, 4, 703, 15, 696, 16, 705, 27, 703, 8, 696, 0, 688, 17, 685, 38, 709,
    0, 674, 7, 698, 5, 694, 18, 694, 0, 673,
];
const vBLS: u32 = 661;
const matrixBLS: [[u32; 3]; 3] = [[2, 3, 1], [1, 2, 3], [1, 1, 2]]; //AES submatrix
const constantsBLS: [[u32; 3]; 6] = [
    [3u32.pow(100), 2u32.pow(100), 5u32.pow(50)],
    [3u32.pow(110), 2u32.pow(110), 5u32.pow(60)],
    [3u32.pow(120), 2u32.pow(120), 5u32.pow(70)],
    [3u32.pow(130), 2u32.pow(130), 5u32.pow(80)],
    [3u32.pow(140), 2u32.pow(140), 5u32.pow(90)],
    [3u32.pow(150), 2u32.pow(150), 5u32.pow(100)],
];
const SboxBLS: [u32; 661] = [
    248, 131, 29, 335, 367, 614, 598, 627, 611, 505, 478, 601, 294, 251, 37, 518, 241, 303, 486,
    470, 655, 149, 498, 249, 141, 180, 263, 24, 320, 336, 512, 166, 626, 328, 439, 386, 256, 32,
    284, 169, 183, 377, 426, 554, 495, 503, 430, 195, 177, 223, 447, 623, 17, 361, 136, 488, 436,
    621, 88, 425, 62, 9, 261, 485, 395, 27, 204, 266, 334, 213, 127, 528, 456, 135, 143, 270, 159,
    574, 584, 20, 156, 418, 450, 496, 16, 420, 171, 549, 19, 288, 365, 42, 452, 311, 273, 69, 198,
    375, 500, 317, 128, 362, 3, 523, 583, 78, 242, 638, 658, 374, 91, 588, 594, 557, 407, 369, 642,
    304, 593, 605, 298, 499, 575, 221, 567, 134, 422, 534, 276, 85, 644, 48, 353, 95, 97, 71, 77,
    449, 7, 105, 112, 542, 615, 191, 525, 419, 210, 539, 178, 411, 253, 50, 99, 193, 527, 378, 417,
    247, 348, 451, 489, 654, 137, 206, 55, 200, 299, 651, 330, 132, 168, 543, 36, 104, 297, 586,
    455, 116, 145, 120, 519, 578, 170, 535, 564, 465, 129, 568, 467, 103, 218, 329, 444, 121, 281,
    279, 544, 580, 51, 610, 63, 82, 250, 60, 538, 364, 553, 227, 53, 25, 264, 603, 453, 622, 565,
    631, 45, 315, 438, 647, 291, 235, 59, 117, 277, 555, 476, 106, 597, 255, 537, 632, 243, 148,
    641, 464, 569, 548, 415, 175, 309, 272, 228, 211, 41, 345, 283, 501, 526, 524, 620, 589, 108,
    114, 139, 239, 428, 494, 344, 356, 196, 458, 560, 323, 352, 618, 262, 383, 645, 199, 174, 173,
    231, 164, 566, 482, 33, 370, 209, 80, 154, 161, 347, 312, 43, 212, 531, 401, 507, 331, 300,
    639, 205, 595, 236, 636, 333, 338, 454, 44, 86, 271, 305, 462, 475, 562, 533, 619, 520, 529,
    545, 215, 442, 38, 237, 340, 440, 617, 355, 606, 125, 325, 151, 257, 319, 31, 101, 371, 635,
    34, 269, 189, 324, 318, 380, 93, 637, 511, 187, 550, 646, 441, 573, 437, 612, 497, 473, 313,
    346, 460, 219, 372, 163, 252, 181, 282, 293, 26, 625, 376, 138, 546, 342, 468, 506, 278, 186,
    100, 207, 513, 435, 230, 83, 286, 480, 469, 35, 350, 387, 233, 92, 310, 399, 648, 459, 547,
    258, 492, 302, 379, 385, 607, 445, 267, 656, 192, 516, 521, 295, 2, 405, 388, 461, 54, 98, 220,
    244, 245, 587, 46, 111, 260, 416, 49, 214, 591, 4, 130, 570, 240, 203, 592, 75, 343, 424, 289,
    142, 366, 90, 423, 179, 600, 0, 110, 392, 608, 57, 630, 70, 185, 609, 502, 190, 222, 472, 397,
    8, 162, 65, 341, 158, 474, 659, 224, 58, 22, 102, 373, 481, 118, 660, 109, 510, 541, 165, 11,
    479, 226, 280, 650, 74, 412, 484, 434, 351, 363, 194, 146, 556, 457, 184, 559, 409, 1, 28, 477,
    446, 624, 532, 216, 337, 268, 431, 72, 536, 321, 144, 47, 23, 483, 150, 339, 76, 599, 396, 66,
    390, 56, 246, 15, 18, 581, 522, 326, 6, 152, 113, 432, 322, 406, 384, 52, 410, 585, 571, 493,
    68, 254, 39, 596, 602, 463, 653, 13, 327, 628, 613, 21, 582, 634, 115, 358, 515, 572, 393, 359,
    12, 490, 308, 67, 208, 398, 201, 176, 122, 403, 316, 558, 155, 14, 487, 332, 590, 197, 259,
    427, 153, 466, 73, 119, 94, 225, 640, 234, 64, 265, 123, 306, 160, 354, 633, 84, 301, 87, 357,
    389, 429, 652, 577, 287, 561, 404, 229, 307, 182, 402, 290, 381, 126, 274, 540, 238, 563, 172,
    514, 147, 448, 275, 517, 81, 368, 530, 296, 643, 508, 579, 167, 124, 96, 433, 107, 89, 157,
    400, 391, 414, 649, 79, 40, 552, 382, 421, 133, 604, 30, 232, 443, 188, 629, 314, 349, 5, 551,
    360, 394, 10, 217, 509, 657, 285, 202, 292, 140, 471, 491, 413, 408, 616, 61, 576, 504,
];

// Convert representation from tupple in (Z_{s1} x ... x Z_{s_n}) to single element
fn compute_whole_representation(decomposition: [u32; 54]) -> BigUint {
    let mut single: BigUint = num_bigint::BigUint::new(vec![0]);
    // Note that decomposition[53] is s_1, not s_n, so decomposition[2] is s_n
    for k in 27..0 {
        single += decomposition[2*k];
        single *= decomposition[2*k+1];
    }
    single
}

fn small_s_box(x: u32, v: u32) -> u32 {
    let mut y = 0;
    if v == vBLS {
        if x < v {y = SboxBLS[x as usize]} else {y = x}
    } else {
        if x < v {y = x.pow(v - 2) % v} else {y = x}
    }
    y
}

// Lookup-table-based Sbox
fn bar(state: [u32; 3], decomposition: [u32; 54], v: u32) -> [u32; 3] {
    let mut new_state = [0; 3];
    let mut nibbles = [0; 27];
    for i in 0..3 {
        // 1. Decomposition
        // get state value that we are decomposing and manipulating
        let mut intermediate = state[i];
        for k in 0..27 {
            let value = intermediate % decomposition[2*k+1];
            // reduce intermediate representation
            if k < 26 {
                intermediate = (intermediate - value)/decomposition[2*k+3];
            }
            // 2. S-box
            nibbles[k] = small_s_box(value, v);
        }

        // 3. Composition
        for k in 26..=0 {
            new_state[i] += nibbles[k];
            new_state[i] *= decomposition[2*k+1];
        }
    }
    new_state
}

// Element-wise power function
fn brick(state: [u32; 3], order: u32) -> [u32; 3] {
    let mut new_state = [0, 0, 0];
    let z_squared = state[2].pow(2);
    new_state[0] = (z_squared.pow(2)*state[2]) % order;
    new_state[1] = state[0] * (z_squared + state[2] + 2) % order;
    new_state[2] = state[1] * (state[0] * state[0] + 3 * state[0] + 4) % order;
    new_state
}

// Apply affine transformation to state
fn concrete(state: [u32; 3], order: u32, matrix: [[u32; 3]; 3], constants: [u32; 3]) -> [u32; 3] {
    let mut new_state = constants;
    // matrix multiplication
    for i in 0..3 {
        for j in 0..3 {
            new_state[i] += matrix[i][j] * state[j] % order;
        }
    }
    new_state
}

fn full_hash(
    state: [u32; 3],
    decomposition: [u32; 54],
    v: u32,
    matrix: [[u32; 3]; 3],
    constants: [[u32; 3]; 6],
) -> [u32; 3] {
    let order = compute_whole_representation(decomposition);
    let mut new_state = concrete(state, order, matrix, constants[0]);
    new_state = brick(new_state, order);
    new_state = concrete(new_state, order, matrix, constants[1]);
    new_state = brick(new_state, order);
    new_state = concrete(new_state, order, matrix, constants[2]);
    new_state = bar(new_state, decomposition, v);
    new_state = concrete(new_state, order, matrix, constants[3]);
    new_state = brick(new_state, order);
    new_state = concrete(new_state, order, matrix, constants[4]);
    new_state = brick(new_state, order);
    new_state = concrete(new_state, order, matrix, constants[5]);
    new_state
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test1() {
        let input1 = [compute_whole_representation(decompositionBLS)-1; 3];
        let output1 = bar(input1, decompositionBLS, vBLS);
        assert_eq!(input1, output1);
    }

    #[test]
    fn test2() {
        let input2 = [compute_whole_representation(decompositionBLS)-2; 3];
        let output2 = bar(input2, decompositionBLS, vBLS);
        assert_eq!(input2, output2);
    }

    #[test]
    fn test3() {
        let input3 = [0; 3];
        let output3 = bar(input3, decompositionBLS, vBLS);
        assert_eq!(input3, output3);
    }

    #[test]
    fn test4() {
        let input4 = [1; 3];
        let output4 = bar(input4, decompositionBLS, vBLS);
        assert_eq!(input4, output4);
    }
    
    #[test]
    fn test5() {
        let input5 = [5; 3];
        let output5 = bar(input5, decompositionBLS, vBLS);
        assert_eq!(input5, output5);
    }
}

/*
#not working with current sboxes
def test(order,decomposition,v,matrix,constants):

BLS_order = compute_whole_representation(decompositionBLS);
print("BLS: ",BLS_order);

#test(BLS_order,decompositionBLS,vBLS,matrixBLS,constantsBLS);
input = [1,2,3]
print(input)
output=fullHash(input,decompositionBLS,vBLS,matrixBLS,constantsBLS)
print(output);
*/
