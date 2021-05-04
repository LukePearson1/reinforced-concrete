// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Here the constants and fundamental building blocks, which make up the moving
//! parts of the hash function Reinforced Concrete, are defined.

extern crate dusk_bls12_381 as BLS;
use bigint::U256 as u256;
use BLS::BlsScalar as Scalar;

pub const vBLS: Scalar = Scalar([661, 0, 0, 0]);
pub const vu256: u256 = u256([661, 0, 0, 0]);
// Elements of the MDS matrix used; this is the matrix [[2,1,1],[1,2,1],[1,1,2]]
// converted to BLS Scalar form (i.e. in Montgomery form)
pub const matrixBLS: [[Scalar; 3]; 3] = [
    [
        Scalar([
            17179869180,
            12756850513266774020,
            3681868479150465002,
            3479420709561305823,
        ]),
        Scalar([
            8589934590,
            6378425256633387010,
            11064306276430008309,
            1739710354780652911,
        ]),
        Scalar([
            8589934590,
            6378425256633387010,
            11064306276430008309,
            1739710354780652911,
        ]),
    ],
    [
        Scalar([
            8589934590,
            6378425256633387010,
            11064306276430008309,
            1739710354780652911,
        ]),
        Scalar([
            17179869180,
            12756850513266774020,
            3681868479150465002,
            3479420709561305823,
        ]),
        Scalar([
            8589934590,
            6378425256633387010,
            11064306276430008309,
            1739710354780652911,
        ]),
    ],
    [
        Scalar([
            8589934590,
            6378425256633387010,
            11064306276430008309,
            1739710354780652911,
        ]),
        Scalar([
            8589934590,
            6378425256633387010,
            11064306276430008309,
            1739710354780652911,
        ]),
        Scalar([
            17179869180,
            12756850513266774020,
            3681868479150465002,
            3479420709561305823,
        ]),
    ],
];

// Constant round vector that is included in concrete (in Montgomery form)
pub const constantsBLS: [[Scalar; 3]; 6] = [
    [
        Scalar([
            15103436521110050711,
            7331784790172381998,
            12044656176741403427,
            3751008647294553039,
        ]),
        Scalar([
            11347474211269653648,
            4128803420964107347,
            6214707797160115797,
            4112634405610528689,
        ]),
        Scalar([
            17304014358209917584,
            2208272603974818668,
            10268468444210955228,
            6179877239278854801,
        ]),
    ],
    [
        Scalar([
            16534261342164866333,
            6832202303744756170,
            1684251025634136077,
            8163602255652314233,
        ]),
        Scalar([
            16811572141480885768,
            6086740446641139821,
            2480905392383068776,
            1165134175098904372,
        ]),
        Scalar([
            17107702180119466287,
            4822673964486880457,
            9690930597399869386,
            5205372976176261225,
        ]),
    ],
    [
        Scalar([
            774652113136442555,
            16405992112841771651,
            6017565519437953844,
            4505701757988831749,
        ]),
        Scalar([
            4237652715300724594,
            7969055898639970355,
            5598934159073726045,
            6898001257326269053,
        ]),
        Scalar([
            2313248965582223868,
            13576718295572789794,
            11585600013457612589,
            5850954271136432218,
        ]),
    ],
    [
        Scalar([
            12954206658231352586,
            4453253639687023763,
            9511288471238421433,
            6024650399278547742,
        ]),
        Scalar([
            4371526775444718771,
            17735767830574081080,
            13237695147841645604,
            4831541254639807973,
        ]),
        Scalar([
            14891113154929442632,
            3849496164829441429,
            2848154445163346576,
            2914992698027003878,
        ]),
    ],
    [
        Scalar([
            1812640294639214912,
            12716517994398543712,
            11407848676021953699,
            6712449845925535010,
        ]),
        Scalar([
            12331354760301169072,
            16276070422768526146,
            7003893990016485254,
            2216263948209348584,
        ]),
        Scalar([
            12276154328438218111,
            16836562938552946063,
            6260598915119486614,
            3428279404074293598,
        ]),
    ],
    [
        Scalar([
            6587846275791075816,
            7625938549621093986,
            13183807980954838340,
            5783004187723978816,
        ]),
        Scalar([
            9734329294999961329,
            15789237201575454651,
            10454603626512543203,
            5651214051507175958,
        ]),
        Scalar([
            18253729846641067784,
            12726484722006391827,
            14165217671517419137,
            7711433387056390051,
        ]),
    ],
];

// decomposition = [v_n, v_{n-1} ..., v_1], which is the representation of q-1
pub const BLS_scalar_decomposition: [u256; 27] = [
    u256([32, 0, 0, 0]),
    u256([35, 0, 0, 0]),
    u256([20, 0, 0, 0]),
    u256([4, 0, 0, 0]),
    u256([4, 0, 0, 0]),
    u256([12, 0, 0, 0]),
    u256([3, 0, 0, 0]),
    u256([9, 0, 0, 0]),
    u256([30, 0, 0, 0]),
    u256([3, 0, 0, 0]),
    u256([27, 0, 0, 0]),
    u256([27, 0, 0, 0]),
    u256([20, 0, 0, 0]),
    u256([31, 0, 0, 0]),
    u256([12, 0, 0, 0]),
    u256([32, 0, 0, 0]),
    u256([39, 0, 0, 0]),
    u256([22, 0, 0, 0]),
    u256([2, 0, 0, 0]),
    u256([11, 0, 0, 0]),
    u256([7, 0, 0, 0]),
    u256([42, 0, 0, 0]),
    u256([11, 0, 0, 0]),
    u256([3, 0, 0, 0]),
    u256([4, 0, 0, 0]),
    u256([28, 0, 0, 0]),
    u256([0, 0, 0, 0]),
];

// decomposition = [s_n, s_{n-1} ..., s_1]
pub const decomposition_s_i: [u256; 27] = [
    u256([693, 0, 0, 0]),
    u256([696, 0, 0, 0]),
    u256([694, 0, 0, 0]),
    u256([668, 0, 0, 0]),
    u256([679, 0, 0, 0]),
    u256([695, 0, 0, 0]),
    u256([691, 0, 0, 0]),
    u256([693, 0, 0, 0]),
    u256([700, 0, 0, 0]),
    u256([688, 0, 0, 0]),
    u256([700, 0, 0, 0]),
    u256([694, 0, 0, 0]),
    u256([701, 0, 0, 0]),
    u256([694, 0, 0, 0]),
    u256([699, 0, 0, 0]),
    u256([701, 0, 0, 0]),
    u256([701, 0, 0, 0]),
    u256([701, 0, 0, 0]),
    u256([695, 0, 0, 0]),
    u256([698, 0, 0, 0]),
    u256([697, 0, 0, 0]),
    u256([703, 0, 0, 0]),
    u256([702, 0, 0, 0]),
    u256([691, 0, 0, 0]),
    u256([688, 0, 0, 0]),
    u256([703, 0, 0, 0]),
    u256([679, 0, 0, 0]),
];

// decomposition_inverses (are not in Montgomery form) = [s_n^{-1}, ...,
// s_1^{-1}]
pub const inverses_s_i: [Scalar; 27] = [
    Scalar([
        2602974132398977206,
        836676680507352349,
        1624520412375634914,
        784276341631287845,
    ]),
    Scalar([
        14439016772699306960,
        4880053042494532754,
        425405000111977086,
        3301353513273012702,
    ]),
    Scalar([
        11184667735571951804,
        2530713830546604934,
        9464048513681697485,
        3792338212704027110,
    ]),
    Scalar([
        17584814563665354469,
        11931036503688296974,
        16968017926336086938,
        3589796777771841665,
    ]),
    Scalar([
        12545619903856777888,
        12888235771846388709,
        10863565764890696465,
        6139808279466941819,
    ]),
    Scalar([
        17512131283132517804,
        5442999654264100671,
        13623050306846729539,
        4700359713576824410,
    ]),
    Scalar([
        6027560511247117796,
        17067218231336313649,
        13728930705354457222,
        6178251843808432543,
    ]),
    Scalar([
        2602974132398977206,
        836676680507352349,
        1624520412375634914,
        784276341631287845,
    ]),
    Scalar([
        2761411830179995378,
        2383420725998368036,
        8729429052584952794,
        3950769984811465720,
    ]),
    Scalar([
        12032948026114369833,
        9483378357172045150,
        13815444268081057867,
        5476687518174069660,
    ]),
    Scalar([
        2761411830179995378,
        2383420725998368036,
        8729429052584952794,
        3950769984811465720,
    ]),
    Scalar([
        11184667735571951804,
        2530713830546604934,
        9464048513681697485,
        3792338212704027110,
    ]),
    Scalar([
        7862562955693360704,
        2189145776794569758,
        17637406801325622816,
        2515145458248633497,
    ]),
    Scalar([
        11184667735571951804,
        2530713830546604934,
        9464048513681697485,
        3792338212704027110,
    ]),
    Scalar([
        14456217320622407667,
        15720950765214895364,
        9290001455579958162,
        239764300535883055,
    ]),
    Scalar([
        7862562955693360704,
        2189145776794569758,
        17637406801325622816,
        2515145458248633497,
    ]),
    Scalar([
        7862562955693360704,
        2189145776794569758,
        17637406801325622816,
        2515145458248633497,
    ]),
    Scalar([
        7862562955693360704,
        2189145776794569758,
        17637406801325622816,
        2515145458248633497,
    ]),
    Scalar([
        17512131283132517804,
        5442999654264100671,
        13623050306846729539,
        4700359713576824410,
    ]),
    Scalar([
        2293620184373320964,
        11321873130585980605,
        6174511658688624013,
        1472790089683983580,
    ]),
    Scalar([
        2667433723633177098,
        9454358102475754882,
        4870880887501134046,
        2457669820768300410,
    ]),
    Scalar([
        18205007880878314118,
        13252090256230760799,
        8933760750817715047,
        4718166574811777612,
    ]),
    Scalar([
        1860112511821192083,
        14052111655569484648,
        3827050090417538583,
        2059377956656186619,
    ]),
    Scalar([
        6027560511247117796,
        17067218231336313649,
        13728930705354457222,
        6178251843808432543,
    ]),
    Scalar([
        12032948026114369833,
        9483378357172045150,
        13815444268081057867,
        5476687518174069660,
    ]),
    Scalar([
        18205007880878314118,
        13252090256230760799,
        8933760750817715047,
        4718166574811777612,
    ]),
    Scalar([
        12545619903856777888,
        12888235771846388709,
        10863565764890696465,
        6139808279466941819,
    ]),
];

pub const SboxBLS: [u256; 661] = [
    u256([248, 0, 0, 0]),
    u256([131, 0, 0, 0]),
    u256([29, 0, 0, 0]),
    u256([335, 0, 0, 0]),
    u256([367, 0, 0, 0]),
    u256([614, 0, 0, 0]),
    u256([598, 0, 0, 0]),
    u256([627, 0, 0, 0]),
    u256([611, 0, 0, 0]),
    u256([505, 0, 0, 0]),
    u256([478, 0, 0, 0]),
    u256([601, 0, 0, 0]),
    u256([294, 0, 0, 0]),
    u256([251, 0, 0, 0]),
    u256([37, 0, 0, 0]),
    u256([518, 0, 0, 0]),
    u256([241, 0, 0, 0]),
    u256([303, 0, 0, 0]),
    u256([486, 0, 0, 0]),
    u256([470, 0, 0, 0]),
    u256([655, 0, 0, 0]),
    u256([149, 0, 0, 0]),
    u256([498, 0, 0, 0]),
    u256([249, 0, 0, 0]),
    u256([141, 0, 0, 0]),
    u256([180, 0, 0, 0]),
    u256([263, 0, 0, 0]),
    u256([24, 0, 0, 0]),
    u256([320, 0, 0, 0]),
    u256([336, 0, 0, 0]),
    u256([512, 0, 0, 0]),
    u256([166, 0, 0, 0]),
    u256([626, 0, 0, 0]),
    u256([328, 0, 0, 0]),
    u256([439, 0, 0, 0]),
    u256([386, 0, 0, 0]),
    u256([256, 0, 0, 0]),
    u256([32, 0, 0, 0]),
    u256([284, 0, 0, 0]),
    u256([169, 0, 0, 0]),
    u256([183, 0, 0, 0]),
    u256([377, 0, 0, 0]),
    u256([426, 0, 0, 0]),
    u256([554, 0, 0, 0]),
    u256([495, 0, 0, 0]),
    u256([503, 0, 0, 0]),
    u256([430, 0, 0, 0]),
    u256([195, 0, 0, 0]),
    u256([177, 0, 0, 0]),
    u256([223, 0, 0, 0]),
    u256([447, 0, 0, 0]),
    u256([623, 0, 0, 0]),
    u256([17, 0, 0, 0]),
    u256([361, 0, 0, 0]),
    u256([136, 0, 0, 0]),
    u256([488, 0, 0, 0]),
    u256([436, 0, 0, 0]),
    u256([621, 0, 0, 0]),
    u256([88, 0, 0, 0]),
    u256([425, 0, 0, 0]),
    u256([62, 0, 0, 0]),
    u256([9, 0, 0, 0]),
    u256([261, 0, 0, 0]),
    u256([485, 0, 0, 0]),
    u256([395, 0, 0, 0]),
    u256([27, 0, 0, 0]),
    u256([204, 0, 0, 0]),
    u256([266, 0, 0, 0]),
    u256([334, 0, 0, 0]),
    u256([213, 0, 0, 0]),
    u256([127, 0, 0, 0]),
    u256([528, 0, 0, 0]),
    u256([456, 0, 0, 0]),
    u256([135, 0, 0, 0]),
    u256([143, 0, 0, 0]),
    u256([270, 0, 0, 0]),
    u256([159, 0, 0, 0]),
    u256([574, 0, 0, 0]),
    u256([584, 0, 0, 0]),
    u256([20, 0, 0, 0]),
    u256([156, 0, 0, 0]),
    u256([418, 0, 0, 0]),
    u256([450, 0, 0, 0]),
    u256([496, 0, 0, 0]),
    u256([16, 0, 0, 0]),
    u256([420, 0, 0, 0]),
    u256([171, 0, 0, 0]),
    u256([549, 0, 0, 0]),
    u256([19, 0, 0, 0]),
    u256([288, 0, 0, 0]),
    u256([365, 0, 0, 0]),
    u256([42, 0, 0, 0]),
    u256([452, 0, 0, 0]),
    u256([311, 0, 0, 0]),
    u256([273, 0, 0, 0]),
    u256([69, 0, 0, 0]),
    u256([198, 0, 0, 0]),
    u256([375, 0, 0, 0]),
    u256([500, 0, 0, 0]),
    u256([317, 0, 0, 0]),
    u256([128, 0, 0, 0]),
    u256([362, 0, 0, 0]),
    u256([3, 0, 0, 0]),
    u256([523, 0, 0, 0]),
    u256([583, 0, 0, 0]),
    u256([78, 0, 0, 0]),
    u256([242, 0, 0, 0]),
    u256([638, 0, 0, 0]),
    u256([658, 0, 0, 0]),
    u256([374, 0, 0, 0]),
    u256([91, 0, 0, 0]),
    u256([588, 0, 0, 0]),
    u256([594, 0, 0, 0]),
    u256([557, 0, 0, 0]),
    u256([407, 0, 0, 0]),
    u256([369, 0, 0, 0]),
    u256([642, 0, 0, 0]),
    u256([304, 0, 0, 0]),
    u256([593, 0, 0, 0]),
    u256([605, 0, 0, 0]),
    u256([298, 0, 0, 0]),
    u256([499, 0, 0, 0]),
    u256([575, 0, 0, 0]),
    u256([221, 0, 0, 0]),
    u256([567, 0, 0, 0]),
    u256([134, 0, 0, 0]),
    u256([422, 0, 0, 0]),
    u256([534, 0, 0, 0]),
    u256([276, 0, 0, 0]),
    u256([85, 0, 0, 0]),
    u256([644, 0, 0, 0]),
    u256([48, 0, 0, 0]),
    u256([353, 0, 0, 0]),
    u256([95, 0, 0, 0]),
    u256([97, 0, 0, 0]),
    u256([71, 0, 0, 0]),
    u256([77, 0, 0, 0]),
    u256([449, 0, 0, 0]),
    u256([7, 0, 0, 0]),
    u256([105, 0, 0, 0]),
    u256([112, 0, 0, 0]),
    u256([542, 0, 0, 0]),
    u256([615, 0, 0, 0]),
    u256([191, 0, 0, 0]),
    u256([525, 0, 0, 0]),
    u256([419, 0, 0, 0]),
    u256([210, 0, 0, 0]),
    u256([539, 0, 0, 0]),
    u256([178, 0, 0, 0]),
    u256([411, 0, 0, 0]),
    u256([253, 0, 0, 0]),
    u256([50, 0, 0, 0]),
    u256([99, 0, 0, 0]),
    u256([193, 0, 0, 0]),
    u256([527, 0, 0, 0]),
    u256([378, 0, 0, 0]),
    u256([417, 0, 0, 0]),
    u256([247, 0, 0, 0]),
    u256([348, 0, 0, 0]),
    u256([451, 0, 0, 0]),
    u256([489, 0, 0, 0]),
    u256([654, 0, 0, 0]),
    u256([137, 0, 0, 0]),
    u256([206, 0, 0, 0]),
    u256([55, 0, 0, 0]),
    u256([200, 0, 0, 0]),
    u256([299, 0, 0, 0]),
    u256([651, 0, 0, 0]),
    u256([330, 0, 0, 0]),
    u256([132, 0, 0, 0]),
    u256([168, 0, 0, 0]),
    u256([543, 0, 0, 0]),
    u256([36, 0, 0, 0]),
    u256([104, 0, 0, 0]),
    u256([297, 0, 0, 0]),
    u256([586, 0, 0, 0]),
    u256([455, 0, 0, 0]),
    u256([116, 0, 0, 0]),
    u256([145, 0, 0, 0]),
    u256([120, 0, 0, 0]),
    u256([519, 0, 0, 0]),
    u256([578, 0, 0, 0]),
    u256([170, 0, 0, 0]),
    u256([535, 0, 0, 0]),
    u256([564, 0, 0, 0]),
    u256([465, 0, 0, 0]),
    u256([129, 0, 0, 0]),
    u256([568, 0, 0, 0]),
    u256([467, 0, 0, 0]),
    u256([103, 0, 0, 0]),
    u256([218, 0, 0, 0]),
    u256([329, 0, 0, 0]),
    u256([444, 0, 0, 0]),
    u256([121, 0, 0, 0]),
    u256([281, 0, 0, 0]),
    u256([279, 0, 0, 0]),
    u256([544, 0, 0, 0]),
    u256([580, 0, 0, 0]),
    u256([51, 0, 0, 0]),
    u256([610, 0, 0, 0]),
    u256([63, 0, 0, 0]),
    u256([82, 0, 0, 0]),
    u256([250, 0, 0, 0]),
    u256([60, 0, 0, 0]),
    u256([538, 0, 0, 0]),
    u256([364, 0, 0, 0]),
    u256([553, 0, 0, 0]),
    u256([227, 0, 0, 0]),
    u256([53, 0, 0, 0]),
    u256([25, 0, 0, 0]),
    u256([264, 0, 0, 0]),
    u256([603, 0, 0, 0]),
    u256([453, 0, 0, 0]),
    u256([622, 0, 0, 0]),
    u256([565, 0, 0, 0]),
    u256([631, 0, 0, 0]),
    u256([45, 0, 0, 0]),
    u256([315, 0, 0, 0]),
    u256([438, 0, 0, 0]),
    u256([647, 0, 0, 0]),
    u256([291, 0, 0, 0]),
    u256([235, 0, 0, 0]),
    u256([59, 0, 0, 0]),
    u256([117, 0, 0, 0]),
    u256([277, 0, 0, 0]),
    u256([555, 0, 0, 0]),
    u256([476, 0, 0, 0]),
    u256([106, 0, 0, 0]),
    u256([597, 0, 0, 0]),
    u256([255, 0, 0, 0]),
    u256([537, 0, 0, 0]),
    u256([632, 0, 0, 0]),
    u256([243, 0, 0, 0]),
    u256([148, 0, 0, 0]),
    u256([641, 0, 0, 0]),
    u256([464, 0, 0, 0]),
    u256([569, 0, 0, 0]),
    u256([548, 0, 0, 0]),
    u256([415, 0, 0, 0]),
    u256([175, 0, 0, 0]),
    u256([309, 0, 0, 0]),
    u256([272, 0, 0, 0]),
    u256([228, 0, 0, 0]),
    u256([211, 0, 0, 0]),
    u256([41, 0, 0, 0]),
    u256([345, 0, 0, 0]),
    u256([283, 0, 0, 0]),
    u256([501, 0, 0, 0]),
    u256([526, 0, 0, 0]),
    u256([524, 0, 0, 0]),
    u256([620, 0, 0, 0]),
    u256([589, 0, 0, 0]),
    u256([108, 0, 0, 0]),
    u256([114, 0, 0, 0]),
    u256([139, 0, 0, 0]),
    u256([239, 0, 0, 0]),
    u256([428, 0, 0, 0]),
    u256([494, 0, 0, 0]),
    u256([344, 0, 0, 0]),
    u256([356, 0, 0, 0]),
    u256([196, 0, 0, 0]),
    u256([458, 0, 0, 0]),
    u256([560, 0, 0, 0]),
    u256([323, 0, 0, 0]),
    u256([352, 0, 0, 0]),
    u256([618, 0, 0, 0]),
    u256([262, 0, 0, 0]),
    u256([383, 0, 0, 0]),
    u256([645, 0, 0, 0]),
    u256([199, 0, 0, 0]),
    u256([174, 0, 0, 0]),
    u256([173, 0, 0, 0]),
    u256([231, 0, 0, 0]),
    u256([164, 0, 0, 0]),
    u256([566, 0, 0, 0]),
    u256([482, 0, 0, 0]),
    u256([33, 0, 0, 0]),
    u256([370, 0, 0, 0]),
    u256([209, 0, 0, 0]),
    u256([80, 0, 0, 0]),
    u256([154, 0, 0, 0]),
    u256([161, 0, 0, 0]),
    u256([347, 0, 0, 0]),
    u256([312, 0, 0, 0]),
    u256([43, 0, 0, 0]),
    u256([212, 0, 0, 0]),
    u256([531, 0, 0, 0]),
    u256([401, 0, 0, 0]),
    u256([507, 0, 0, 0]),
    u256([331, 0, 0, 0]),
    u256([300, 0, 0, 0]),
    u256([639, 0, 0, 0]),
    u256([205, 0, 0, 0]),
    u256([595, 0, 0, 0]),
    u256([236, 0, 0, 0]),
    u256([636, 0, 0, 0]),
    u256([333, 0, 0, 0]),
    u256([338, 0, 0, 0]),
    u256([454, 0, 0, 0]),
    u256([44, 0, 0, 0]),
    u256([86, 0, 0, 0]),
    u256([271, 0, 0, 0]),
    u256([305, 0, 0, 0]),
    u256([462, 0, 0, 0]),
    u256([475, 0, 0, 0]),
    u256([562, 0, 0, 0]),
    u256([533, 0, 0, 0]),
    u256([619, 0, 0, 0]),
    u256([520, 0, 0, 0]),
    u256([529, 0, 0, 0]),
    u256([545, 0, 0, 0]),
    u256([215, 0, 0, 0]),
    u256([442, 0, 0, 0]),
    u256([38, 0, 0, 0]),
    u256([237, 0, 0, 0]),
    u256([340, 0, 0, 0]),
    u256([440, 0, 0, 0]),
    u256([617, 0, 0, 0]),
    u256([355, 0, 0, 0]),
    u256([606, 0, 0, 0]),
    u256([125, 0, 0, 0]),
    u256([325, 0, 0, 0]),
    u256([151, 0, 0, 0]),
    u256([257, 0, 0, 0]),
    u256([319, 0, 0, 0]),
    u256([31, 0, 0, 0]),
    u256([101, 0, 0, 0]),
    u256([371, 0, 0, 0]),
    u256([635, 0, 0, 0]),
    u256([34, 0, 0, 0]),
    u256([269, 0, 0, 0]),
    u256([189, 0, 0, 0]),
    u256([324, 0, 0, 0]),
    u256([318, 0, 0, 0]),
    u256([380, 0, 0, 0]),
    u256([93, 0, 0, 0]),
    u256([637, 0, 0, 0]),
    u256([511, 0, 0, 0]),
    u256([187, 0, 0, 0]),
    u256([550, 0, 0, 0]),
    u256([646, 0, 0, 0]),
    u256([441, 0, 0, 0]),
    u256([573, 0, 0, 0]),
    u256([437, 0, 0, 0]),
    u256([612, 0, 0, 0]),
    u256([497, 0, 0, 0]),
    u256([473, 0, 0, 0]),
    u256([313, 0, 0, 0]),
    u256([346, 0, 0, 0]),
    u256([460, 0, 0, 0]),
    u256([219, 0, 0, 0]),
    u256([372, 0, 0, 0]),
    u256([163, 0, 0, 0]),
    u256([252, 0, 0, 0]),
    u256([181, 0, 0, 0]),
    u256([282, 0, 0, 0]),
    u256([293, 0, 0, 0]),
    u256([26, 0, 0, 0]),
    u256([625, 0, 0, 0]),
    u256([376, 0, 0, 0]),
    u256([138, 0, 0, 0]),
    u256([546, 0, 0, 0]),
    u256([342, 0, 0, 0]),
    u256([468, 0, 0, 0]),
    u256([506, 0, 0, 0]),
    u256([278, 0, 0, 0]),
    u256([186, 0, 0, 0]),
    u256([100, 0, 0, 0]),
    u256([207, 0, 0, 0]),
    u256([513, 0, 0, 0]),
    u256([435, 0, 0, 0]),
    u256([230, 0, 0, 0]),
    u256([83, 0, 0, 0]),
    u256([286, 0, 0, 0]),
    u256([480, 0, 0, 0]),
    u256([469, 0, 0, 0]),
    u256([35, 0, 0, 0]),
    u256([350, 0, 0, 0]),
    u256([387, 0, 0, 0]),
    u256([233, 0, 0, 0]),
    u256([92, 0, 0, 0]),
    u256([310, 0, 0, 0]),
    u256([399, 0, 0, 0]),
    u256([648, 0, 0, 0]),
    u256([459, 0, 0, 0]),
    u256([547, 0, 0, 0]),
    u256([258, 0, 0, 0]),
    u256([492, 0, 0, 0]),
    u256([302, 0, 0, 0]),
    u256([379, 0, 0, 0]),
    u256([385, 0, 0, 0]),
    u256([607, 0, 0, 0]),
    u256([445, 0, 0, 0]),
    u256([267, 0, 0, 0]),
    u256([656, 0, 0, 0]),
    u256([192, 0, 0, 0]),
    u256([516, 0, 0, 0]),
    u256([521, 0, 0, 0]),
    u256([295, 0, 0, 0]),
    u256([2, 0, 0, 0]),
    u256([405, 0, 0, 0]),
    u256([388, 0, 0, 0]),
    u256([461, 0, 0, 0]),
    u256([54, 0, 0, 0]),
    u256([98, 0, 0, 0]),
    u256([220, 0, 0, 0]),
    u256([244, 0, 0, 0]),
    u256([245, 0, 0, 0]),
    u256([587, 0, 0, 0]),
    u256([46, 0, 0, 0]),
    u256([111, 0, 0, 0]),
    u256([260, 0, 0, 0]),
    u256([416, 0, 0, 0]),
    u256([49, 0, 0, 0]),
    u256([214, 0, 0, 0]),
    u256([591, 0, 0, 0]),
    u256([4, 0, 0, 0]),
    u256([130, 0, 0, 0]),
    u256([570, 0, 0, 0]),
    u256([240, 0, 0, 0]),
    u256([203, 0, 0, 0]),
    u256([592, 0, 0, 0]),
    u256([75, 0, 0, 0]),
    u256([343, 0, 0, 0]),
    u256([424, 0, 0, 0]),
    u256([289, 0, 0, 0]),
    u256([142, 0, 0, 0]),
    u256([366, 0, 0, 0]),
    u256([90, 0, 0, 0]),
    u256([423, 0, 0, 0]),
    u256([179, 0, 0, 0]),
    u256([600, 0, 0, 0]),
    u256([0, 0, 0, 0]),
    u256([110, 0, 0, 0]),
    u256([392, 0, 0, 0]),
    u256([608, 0, 0, 0]),
    u256([57, 0, 0, 0]),
    u256([630, 0, 0, 0]),
    u256([70, 0, 0, 0]),
    u256([185, 0, 0, 0]),
    u256([609, 0, 0, 0]),
    u256([502, 0, 0, 0]),
    u256([190, 0, 0, 0]),
    u256([222, 0, 0, 0]),
    u256([472, 0, 0, 0]),
    u256([397, 0, 0, 0]),
    u256([8, 0, 0, 0]),
    u256([162, 0, 0, 0]),
    u256([65, 0, 0, 0]),
    u256([341, 0, 0, 0]),
    u256([158, 0, 0, 0]),
    u256([474, 0, 0, 0]),
    u256([659, 0, 0, 0]),
    u256([224, 0, 0, 0]),
    u256([58, 0, 0, 0]),
    u256([22, 0, 0, 0]),
    u256([102, 0, 0, 0]),
    u256([373, 0, 0, 0]),
    u256([481, 0, 0, 0]),
    u256([118, 0, 0, 0]),
    u256([660, 0, 0, 0]),
    u256([109, 0, 0, 0]),
    u256([510, 0, 0, 0]),
    u256([541, 0, 0, 0]),
    u256([165, 0, 0, 0]),
    u256([11, 0, 0, 0]),
    u256([479, 0, 0, 0]),
    u256([226, 0, 0, 0]),
    u256([280, 0, 0, 0]),
    u256([650, 0, 0, 0]),
    u256([74, 0, 0, 0]),
    u256([412, 0, 0, 0]),
    u256([484, 0, 0, 0]),
    u256([434, 0, 0, 0]),
    u256([351, 0, 0, 0]),
    u256([363, 0, 0, 0]),
    u256([194, 0, 0, 0]),
    u256([146, 0, 0, 0]),
    u256([556, 0, 0, 0]),
    u256([457, 0, 0, 0]),
    u256([184, 0, 0, 0]),
    u256([559, 0, 0, 0]),
    u256([409, 0, 0, 0]),
    u256([1, 0, 0, 0]),
    u256([28, 0, 0, 0]),
    u256([477, 0, 0, 0]),
    u256([446, 0, 0, 0]),
    u256([624, 0, 0, 0]),
    u256([532, 0, 0, 0]),
    u256([216, 0, 0, 0]),
    u256([337, 0, 0, 0]),
    u256([268, 0, 0, 0]),
    u256([431, 0, 0, 0]),
    u256([72, 0, 0, 0]),
    u256([536, 0, 0, 0]),
    u256([321, 0, 0, 0]),
    u256([144, 0, 0, 0]),
    u256([47, 0, 0, 0]),
    u256([23, 0, 0, 0]),
    u256([483, 0, 0, 0]),
    u256([150, 0, 0, 0]),
    u256([339, 0, 0, 0]),
    u256([76, 0, 0, 0]),
    u256([599, 0, 0, 0]),
    u256([396, 0, 0, 0]),
    u256([66, 0, 0, 0]),
    u256([390, 0, 0, 0]),
    u256([56, 0, 0, 0]),
    u256([246, 0, 0, 0]),
    u256([15, 0, 0, 0]),
    u256([18, 0, 0, 0]),
    u256([581, 0, 0, 0]),
    u256([522, 0, 0, 0]),
    u256([326, 0, 0, 0]),
    u256([6, 0, 0, 0]),
    u256([152, 0, 0, 0]),
    u256([113, 0, 0, 0]),
    u256([432, 0, 0, 0]),
    u256([322, 0, 0, 0]),
    u256([406, 0, 0, 0]),
    u256([384, 0, 0, 0]),
    u256([52, 0, 0, 0]),
    u256([410, 0, 0, 0]),
    u256([585, 0, 0, 0]),
    u256([571, 0, 0, 0]),
    u256([493, 0, 0, 0]),
    u256([68, 0, 0, 0]),
    u256([254, 0, 0, 0]),
    u256([39, 0, 0, 0]),
    u256([596, 0, 0, 0]),
    u256([602, 0, 0, 0]),
    u256([463, 0, 0, 0]),
    u256([653, 0, 0, 0]),
    u256([13, 0, 0, 0]),
    u256([327, 0, 0, 0]),
    u256([628, 0, 0, 0]),
    u256([613, 0, 0, 0]),
    u256([21, 0, 0, 0]),
    u256([582, 0, 0, 0]),
    u256([634, 0, 0, 0]),
    u256([115, 0, 0, 0]),
    u256([358, 0, 0, 0]),
    u256([515, 0, 0, 0]),
    u256([572, 0, 0, 0]),
    u256([393, 0, 0, 0]),
    u256([359, 0, 0, 0]),
    u256([12, 0, 0, 0]),
    u256([490, 0, 0, 0]),
    u256([308, 0, 0, 0]),
    u256([67, 0, 0, 0]),
    u256([208, 0, 0, 0]),
    u256([398, 0, 0, 0]),
    u256([201, 0, 0, 0]),
    u256([176, 0, 0, 0]),
    u256([122, 0, 0, 0]),
    u256([403, 0, 0, 0]),
    u256([316, 0, 0, 0]),
    u256([558, 0, 0, 0]),
    u256([155, 0, 0, 0]),
    u256([14, 0, 0, 0]),
    u256([487, 0, 0, 0]),
    u256([332, 0, 0, 0]),
    u256([590, 0, 0, 0]),
    u256([197, 0, 0, 0]),
    u256([259, 0, 0, 0]),
    u256([427, 0, 0, 0]),
    u256([153, 0, 0, 0]),
    u256([466, 0, 0, 0]),
    u256([73, 0, 0, 0]),
    u256([119, 0, 0, 0]),
    u256([94, 0, 0, 0]),
    u256([225, 0, 0, 0]),
    u256([640, 0, 0, 0]),
    u256([234, 0, 0, 0]),
    u256([64, 0, 0, 0]),
    u256([265, 0, 0, 0]),
    u256([123, 0, 0, 0]),
    u256([306, 0, 0, 0]),
    u256([160, 0, 0, 0]),
    u256([354, 0, 0, 0]),
    u256([633, 0, 0, 0]),
    u256([84, 0, 0, 0]),
    u256([301, 0, 0, 0]),
    u256([87, 0, 0, 0]),
    u256([357, 0, 0, 0]),
    u256([389, 0, 0, 0]),
    u256([429, 0, 0, 0]),
    u256([652, 0, 0, 0]),
    u256([577, 0, 0, 0]),
    u256([287, 0, 0, 0]),
    u256([561, 0, 0, 0]),
    u256([404, 0, 0, 0]),
    u256([229, 0, 0, 0]),
    u256([307, 0, 0, 0]),
    u256([182, 0, 0, 0]),
    u256([402, 0, 0, 0]),
    u256([290, 0, 0, 0]),
    u256([381, 0, 0, 0]),
    u256([126, 0, 0, 0]),
    u256([274, 0, 0, 0]),
    u256([540, 0, 0, 0]),
    u256([238, 0, 0, 0]),
    u256([563, 0, 0, 0]),
    u256([172, 0, 0, 0]),
    u256([514, 0, 0, 0]),
    u256([147, 0, 0, 0]),
    u256([448, 0, 0, 0]),
    u256([275, 0, 0, 0]),
    u256([517, 0, 0, 0]),
    u256([81, 0, 0, 0]),
    u256([368, 0, 0, 0]),
    u256([530, 0, 0, 0]),
    u256([296, 0, 0, 0]),
    u256([643, 0, 0, 0]),
    u256([508, 0, 0, 0]),
    u256([579, 0, 0, 0]),
    u256([167, 0, 0, 0]),
    u256([124, 0, 0, 0]),
    u256([96, 0, 0, 0]),
    u256([433, 0, 0, 0]),
    u256([107, 0, 0, 0]),
    u256([89, 0, 0, 0]),
    u256([157, 0, 0, 0]),
    u256([400, 0, 0, 0]),
    u256([391, 0, 0, 0]),
    u256([414, 0, 0, 0]),
    u256([649, 0, 0, 0]),
    u256([79, 0, 0, 0]),
    u256([40, 0, 0, 0]),
    u256([552, 0, 0, 0]),
    u256([382, 0, 0, 0]),
    u256([421, 0, 0, 0]),
    u256([133, 0, 0, 0]),
    u256([604, 0, 0, 0]),
    u256([30, 0, 0, 0]),
    u256([232, 0, 0, 0]),
    u256([443, 0, 0, 0]),
    u256([188, 0, 0, 0]),
    u256([629, 0, 0, 0]),
    u256([314, 0, 0, 0]),
    u256([349, 0, 0, 0]),
    u256([5, 0, 0, 0]),
    u256([551, 0, 0, 0]),
    u256([360, 0, 0, 0]),
    u256([394, 0, 0, 0]),
    u256([10, 0, 0, 0]),
    u256([217, 0, 0, 0]),
    u256([509, 0, 0, 0]),
    u256([657, 0, 0, 0]),
    u256([285, 0, 0, 0]),
    u256([202, 0, 0, 0]),
    u256([292, 0, 0, 0]),
    u256([140, 0, 0, 0]),
    u256([471, 0, 0, 0]),
    u256([491, 0, 0, 0]),
    u256([413, 0, 0, 0]),
    u256([408, 0, 0, 0]),
    u256([616, 0, 0, 0]),
    u256([61, 0, 0, 0]),
    u256([576, 0, 0, 0]),
    u256([504, 0, 0, 0]),
];
