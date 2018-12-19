// numeri pseudo-casuali
// Riccardo Bertossa, 01/11/2013
// compilare con gcc -c rnd.c -o rnd.o
#include <stdint.h> // definisce uint32_t e uint64_t
#include <stdio.h> // funzioni input/output
#include <math.h> // funzioni matematiche
#include <stdlib.h>
#include <stdbool.h>

#ifndef ANALISI
#include <gmp.h>
#include <mpfr.h> // per il calcolo di gamma
#include <mpf2mpfr.h>
#endif

// generatore migliore di lcg, usato per inizializzare cmwc4096
#define SHR3 (jz=jsr, jsr^=(jsr<<13), jsr^=(jsr>>17), jsr^=(jsr<<5),jz+jsr)

#define UNI (rnd_single()-0.5)*2
// array di stato dei generatori, globale
static uint32_t Q_[4096], jz,jsr=123;

void set_SHR3_jsr(uint32_t j)  {
    jsr=j;
}

uint32_t rnd_shr3(){
    return SHR3;
}

// tabelle per il generatore (veloce) di numeri con distribuzione gaussiana
#define PI 3.1415926535897932384626433832795
#define LEVELS 128
#define PARAM_R 3.44261985589665212142432

#define PRECISIONE 256

static float nc2[2] = {0.5, 0.5};
static float nc3[3] = {1.0/3.0, 4.0/3.0, 1.0/3.0};
static float nc4[4] = {3.0/8.0, 9.0/8.0, 9.0/8.0, 3.0/8.0};
static float nc5[5] = {14.0/45.0, 64.0/45.0, 8.0/15.0, 64.0/45.0, 14.0/45.0};
static float nc6[6] = {95.0/288.0, 125.0/96.0, 125.0/144.0, 125.0/144.0, 125.0/96.0, 95.0/288.0};

#ifndef ANALISI
static float chi2_040q[201]={
    0,
    0,
    1.023,
    1.871,
    2.754,
    3.656,
    4.571,
    5.495,
    6.423,
    7.358,
    8.297,
    9.239,
    10.182,
    11.130,
    12.080,
    13.031,
    13.983,
    14.939,
    15.894,
    16.851,
    17.810,
    18.770,
    19.730,
    20.691,
    21.653,
    22.616,
    23.580,
    24.545,
    25.509,
    26.475,
    27.443,
    28.409,
    29.376,
    30.345,
    31.314,
    32.283,
    33.252,
    34.223,
    35.193,
    36.164,
    37.134,
    38.106,
    39.078,
    40.050,
    41.022,
    41.996,
    42.969,
    43.943,
    44.916,
    45.890,
    46.865,
    47.838,
    48.813,
    49.788,
    50.765,
    51.740,
    52.715,
    53.692,
    54.668,
    55.645,
    56.621,
    57.598,
    58.574,
    59.552,
    60.529,
    61.507,
    62.485,
    63.463,
    64.441,
    65.419,
    66.397,
    67.375,
    68.354,
    69.332,
    70.312,
    71.291,
    72.271,
    73.250,
    74.230,
    75.209,
    76.189,
    77.168,
    78.149,
    79.129,
    80.110,
    81.089,
    82.070,
    83.051,
    84.031,
    85.012,
    85.993,
    86.974,
    87.955,
    88.937,
    89.918,
    90.899,
    91.882,
    92.863,
    93.845,
    94.826,
    95.809,
    96.790,
    97.772,
    98.755,
    99.737,
    100.720,
    101.703,
    102.685,
    103.668,
    104.650,
    105.633,
    106.617,
    107.599,
    108.582,
    109.566,
    110.548,
    111.532,
    112.515,
    113.499,
    114.481,
    115.465,
    116.449,
    117.433,
    118.416,
    119.400,
    120.384,
    121.368,
    122.352,
    123.336,
    124.320,
    125.304,
    126.289,
    127.273,
    128.257,
    129.241,
    130.227,
    131.211,
    132.196,
    133.180,
    134.164,
    135.150,
    136.135,
    137.119,
    138.105,
    139.089,
    140.074,
    141.060,
    142.045,
    143.029,
    144.015,
    145.000,
    145.986,
    146.971,
    147.957,
    148.942,
    149.928,
    150.914,
    151.899,
    152.885,
    153.870,
    154.856,
    155.843,
    156.828,
    157.814,
    158.799,
    159.786,
    160.772,
    161.759,
    162.744,
    163.730,
    164.717,
    165.702,
    166.689,
    167.675,
    168.662,
    169.649,
    170.634,
    171.621,
    172.608,
    173.594,
    174.581,
    175.568,
    176.553,
    177.540,
    178.527,
    179.514,
    180.501,
    181.488,
    182.475,
    183.462,
    184.449,
    185.436,
    186.423,
    187.410,
    188.397,
    189.384,
    190.371,
    191.358,
    192.345,
    193.332,
    194.319,
};

static float chi2_050q[201]={
    0,
    0,
    1.388,
    2.367,
    3.357,
    4.352,
    5.349,
    6.347,
    7.346,
    8.343,
    9.342,
    10.341,
    11.342,
    12.341,
    13.340,
    14.340,
    15.339,
    16.340,
    17.339,
    18.338,
    19.338,
    20.339,
    21.338,
    22.338,
    23.337,
    24.338,
    25.337,
    26.337,
    27.336,
    28.337,
    29.337,
    30.336,
    31.337,
    32.336,
    33.336,
    34.337,
    35.336,
    36.336,
    37.337,
    38.336,
    39.336,
    40.335,
    41.336,
    42.336,
    43.335,
    44.336,
    45.336,
    46.335,
    47.336,
    48.335,
    49.335,
    50.336,
    51.335,
    52.336,
    53.336,
    54.335,
    55.336,
    56.336,
    57.335,
    58.336,
    59.336,
    60.335,
    61.336,
    62.335,
    63.335,
    64.336,
    65.335,
    66.335,
    67.336,
    68.335,
    69.335,
    70.336,
    71.335,
    72.335,
    73.336,
    74.335,
    75.335,
    76.336,
    77.335,
    78.335,
    79.336,
    80.335,
    81.335,
    82.334,
    83.335,
    84.335,
    85.334,
    86.335,
    87.335,
    88.334,
    89.335,
    90.335,
    91.334,
    92.335,
    93.335,
    94.334,
    95.335,
    96.335,
    97.334,
    98.335,
    99.335,
    100.335,
    101.335,
    102.336,
    103.335,
    104.335,
    105.336,
    106.335,
    107.335,
    108.334,
    109.335,
    110.335,
    111.334,
    112.335,
    113.335,
    114.334,
    115.335,
    116.335,
    117.334,
    118.335,
    119.335,
    120.334,
    121.335,
    122.335,
    123.334,
    124.335,
    125.335,
    126.334,
    127.335,
    128.335,
    129.334,
    130.335,
    131.335,
    132.334,
    133.335,
    134.335,
    135.334,
    136.335,
    137.335,
    138.334,
    139.335,
    140.334,
    141.334,
    142.335,
    143.334,
    144.334,
    145.335,
    146.334,
    147.334,
    148.335,
    149.334,
    150.335,
    151.335,
    152.334,
    153.335,
    154.335,
    155.334,
    156.335,
    157.335,
    158.334,
    159.335,
    160.335,
    161.334,
    162.335,
    163.335,
    164.334,
    165.335,
    166.335,
    167.334,
    168.335,
    169.335,
    170.334,
    171.335,
    172.335,
    173.334,
    174.335,
    175.335,
    176.334,
    177.335,
    178.334,
    179.334,
    180.335,
    181.334,
    182.334,
    183.335,
    184.334,
    185.334,
    186.335,
    187.334,
    188.334,
    189.335,
    190.334,
    191.334,
    192.335,
    193.334,
    194.334,
    195.335,
    196.334,
    197.334,
    198.335,
    199.334,
};



static float chi2_060q[201]={
    0,
    0,
    1.833,
    2.948,
    4.046,
    5.133,
    6.212,
    7.284,
    8.351,
    9.414,
    10.475,
    11.531,
    12.585,
    13.637,
    14.687,
    15.734,
    16.781,
    17.825,
    18.869,
    19.911,
    20.952,
    21.992,
    23.031,
    24.069,
    25.107,
    26.144,
    27.180,
    28.215,
    29.249,
    30.284,
    31.317,
    32.349,
    33.381,
    34.413,
    35.444,
    36.476,
    37.506,
    38.535,
    39.564,
    40.595,
    41.622,
    42.651,
    43.679,
    44.706,
    45.734,
    46.761,
    47.787,
    48.815,
    49.841,
    50.867,
    51.893,
    52.918,
    53.944,
    54.968,
    55.993,
    57.017,
    58.042,
    59.065,
    60.089,
    61.112,
    62.135,
    63.158,
    64.181,
    65.204,
    66.227,
    67.250,
    68.272,
    69.293,
    70.316,
    71.338,
    72.359,
    73.381,
    74.402,
    75.422,
    76.444,
    77.465,
    78.485,
    79.507,
    80.527,
    81.547,
    82.567,
    83.587,
    84.607,
    85.627,
    86.647,
    87.665,
    88.685,
    89.705,
    90.724,
    91.742,
    92.762,
    93.781,
    94.799,
    95.818,
    96.836,
    97.855,
    98.873,
    99.892,
    100.911,
    101.929,
    102.946,
    103.965,
    104.983,
    106.000,
    107.019,
    108.036,
    109.053,
    110.071,
    111.088,
    112.105,
    113.122,
    114.139,
    115.156,
    116.173,
    117.190,
    118.207,
    119.224,
    120.240,
    121.257,
    122.274,
    123.289,
    124.306,
    125.322,
    126.339,
    127.354,
    128.371,
    129.387,
    130.402,
    131.419,
    132.435,
    133.450,
    134.466,
    135.481,
    136.497,
    137.512,
    138.528,
    139.543,
    140.559,
    141.574,
    142.590,
    143.605,
    144.621,
    145.635,
    146.650,
    147.666,
    148.680,
    149.695,
    150.710,
    151.725,
    152.739,
    153.755,
    154.769,
    155.784,
    156.798,
    157.812,
    158.828,
    159.842,
    160.856,
    161.870,
    162.884,
    163.899,
    164.913,
    165.927,
    166.941,
    167.955,
    168.969,
    169.983,
    170.997,
    172.011,
    173.024,
    174.038,
    175.052,
    176.066,
    177.080,
    178.092,
    179.106,
    180.120,
    181.134,
    182.147,
    183.161,
    184.173,
    185.187,
    186.201,
    187.214,
    188.228,
    189.240,
    190.253,
    191.267,
    192.279,
    193.293,
    194.306,
    195.318,
    196.332,
    197.345,
    198.357,
    199.370,
    200.384,
    201.397,
    202.409,
    203.422,
    204.434,
};

static float chi2_095q[201]={
        0,
        3.841,
        5.993,
        7.815,
        9.489,
        11.072,
        12.593,
        14.067,
        15.509,
        16.920,
        18.308,
        19.676,
        21.027,
        22.362,
        23.685,
        24.996,
        26.297,
        27.588,
        28.871,
        30.144,
        31.412,
        32.672,
        33.926,
        35.174,
        36.416,
        37.653,
        38.886,
        40.113,
        41.337,
    42.558,
        43.773,
        44.985,
        46.194,
        47.400,
        48.603,
        49.802,
        50.999,
        52.193,
        53.384,
        54.574,
        55.759,
        56.944,
        58.124,
        59.305,
        60.482,
        61.657,
        62.830,
        64.001,
        65.171,
        66.340,
        67.505,
        68.671,
        69.833,
        70.994,
        72.154,
        73.312,
        74.470,
        75.625,
        76.778,
        77.932,
        79.082,
        80.233,
        81.382,
        82.529,
        83.675,
        84.821,
        85.966,
    87.109,
        88.250,
        89.392,
        90.532,
        91.670,
        92.809,
        93.946,
        95.083,
        96.217,
        97.351,
        98.485,
        99.617,
        100.750,
        101.880,
        103.011,
        104.140,
        105.268,
        106.396,
        107.523,
        108.649,
        109.774,
        110.899,
        112.023,
        113.146,
        114.268,
        115.390,
        116.512,
        117.633,
        118.752,
        119.872,
        120.990,
        122.109,
        123.226,
        124.342,
        125.460,
        126.574,
        127.690,
        128.805,
        129.919,
        131.032,
        132.145,
        133.258,
        134.370,
        135.481,
        136.591,
        137.703,
        138.811,
        139.921,
        141.030,
        142.138,
        143.247,
        144.354,
        145.461,
        146.568,
        147.675,
        148.780,
        149.886,
        150.990,
        152.094,
        153.198,
        154.302,
        155.406,
        156.509,
        157.611,
        158.712,
        159.815,
        160.916,
        162.017,
        163.116,
        164.217,
        165.317,
        166.416,
        167.516,
        168.614,
        169.712,
        170.810,
        171.908,
        173.004,
        174.102,
        175.199,
        176.295,
        177.390,
        178.487,
        179.582,
        180.677,
        181.770,
        182.865,
        183.959,
        185.052,
        186.146,
        187.239,
        188.333,
        189.425,
        190.517,
        191.609,
        192.701,
        193.793,
        194.883,
        195.974,
        197.064,
        198.155,
        199.245,
        200.335,
        201.424,
        202.514,
        203.602,
        204.691,
        205.780,
        206.867,
        207.956,
        209.044,
        210.131,
        211.217,
        212.305,
        213.391,
        214.478,
        215.564,
        216.650,
        217.736,
        218.821,
        219.907,
        220.991,
        222.076,
        223.162,
        224.245,
        225.329,
        226.414,
        227.497,
        228.581,
        229.664,
        230.747,
        231.830,
        232.912,
        233.995,
};
//ifndef ANALISI
#endif

static int tr=0;
static const double ytab[128] = {
  1,
  0.963599693155767599601648,
  0.936282681708371075467596,
  0.913043647992038059994997,
  0.892281650802302723007165,
  0.873243048926853588791902,
  0.855500607885064284178756,
  0.838783605310647223788196,
  0.822907211395262000500971,
  0.807738294696121113398209,
  0.793177011783859210532779,
  0.779146085941703165632235,
  0.765584173909235994800778,
  0.752441559185703801449899,
  0.739677243683338148559379,
  0.727256918354505877929455,
  0.71515150742047704368224,
  0.703336099025817416326938,
  0.691789143446035852786081,
  0.680491841006414333545856,
  0.669427667357706168906295,
  0.658582000058653680158858,
  0.647941821118550808730822,
  0.637495477343144874281302,
  0.627232485257814565783098,
  0.617143370826562488704402,
  0.60721953663260488550992,
  0.597453150951812282171075,
  0.587837054441820525710691,
  0.57836468112670231279185,
  0.569029991074721612247787,
  0.559827412710694817908181,
  0.550751793121055277381678,
  0.541798355031724123112345,
  0.53296265938998758837921,
  0.524240572678992798085186,
  0.515628238249872051959892,
  0.507122051081304589514673,
  0.498718635476584324220121,
  0.490414825289321637409191,
  0.482207646334838690624641,
  0.474094300698249604526618,
  0.466072152694570979243003,
  0.458138716272871964831001,
  0.450291643686926960859244,
  0.442528715280246614831487,
  0.434847830254661896375265,
  0.427246998309562398802494,
  0.419724332054038234673487,
  0.412278040107024620618785,
  0.404906420811488350994666,
  0.397607856498042552075475,
  0.390380808241389489162694,
  0.383223811059883645872926,
  0.376135469514454429473596,
  0.369114453668275139518502,
  0.362159495373033211617374,
  0.355269384851547144351353,
  0.348442967549872469348716,
  0.341679141235013684121299,
  0.33497685331697116146777,
  0.328335098376152396086927,
  0.321752915879208620313063,
  0.315229388068157522215663,
  0.308763638009251922818042,
  0.302354827789479762737332,
  0.296002156849855836591544,
  0.28970486044581049771064,
  0.283462208226012517787131,
  0.277273502921897711534925,
  0.271138079141025273434483,
  0.265055302258161940294889,
  0.259024567398710742630486,
  0.253045298509765859335653,
  0.247116947514696735815539,
  0.24123899354775131610036,
  0.235410942265727656000122,
  0.229632325234302703552581,
  0.22390269938713388746642,
  0.218221646556370595994758,
  0.212588773073736100604209,
  0.207003709441873800636566,
  0.2014661100762032411818,
  0.195975653118110413993186,
  0.190532040320913721116385,
  0.185134997010713432155985,
  0.179784272124962114239656,
  0.174479638332402322954285,
  0.169220892238924751756607,
  0.164007854684927747911202,
  0.158840371140935078126885,
  0.153718312209586570658669,
  0.148641574243696965656861,
  0.143610080091932994685718,
  0.13862377998585103701693,
  0.133682652584647641179354,
  0.128786706197103961086573,
  0.123935980203981742464476,
  0.119130546708718587667411,
  0.114370512449888274521761,
  0.109656021015817761000189,
  0.104987255410354539778514,
  0.100364441029545540132178,
  0.0957878491225781521596647,
  0.0912578008276347102008865,
  0.0867746718955429689981638,
  0.0823388982429574082428027,
  0.0779509825146547141881268,
  0.0736115018847548933890852,
  0.0693211173941802526012133,
  0.0650805852136318737531138,
  0.0608907703485663759720846,
  0.0567526634815385841881903,
  0.0526674019035031697927905,
  0.0486362958602840518780793,
  0.0446608622008724298002803,
  0.0407428680747906046324968,
  0.0368843887869687741283422,
  0.0330878861465051555658772,
  0.0293563174402538296176444,
  0.0256932919361496165715762,
  0.0221033046161115926152897,
  0.0185921027371658126496087,
  0.0151672980106720424679527,
  0.0118394786579823137151507,
  0.00862448441293047096816743,
  0.00554899522081647053916793,
  0.00266962908390250350923154,
};

static const unsigned long ktab[128] = {
  0,
  1611602771,
  1826899878,
  1918584481,
  1969227037,
  2001281515,
  2023368125,
  2039498179,
  2051788381,
  2061460127,
  2069267110,
  2075699398,
  2081089314,
  2085670119,
  2089610331,
  2093034710,
  2096037586,
  2098691595,
  2101053571,
  2103168620,
  2105072996,
  2106796166,
  2108362327,
  2109791536,
  2111100552,
  2112303493,
  2113412330,
  2114437283,
  2115387130,
  2116269447,
  2117090813,
  2117856962,
  2118572919,
  2119243101,
  2119871411,
  2120461303,
  2121015852,
  2121537798,
  2122029592,
  2122493434,
  2122931299,
  2123344971,
  2123736059,
  2124106020,
  2124456175,
  2124787725,
  2125101763,
  2125399283,
  2125681194,
  2125948325,
  2126201433,
  2126441213,
  2126668298,
  2126883268,
  2127086657,
  2127278949,
  2127460589,
  2127631985,
  2127793506,
  2127945490,
  2128088244,
  2128222044,
  2128347141,
  2128463758,
  2128572095,
  2128672327,
  2128764606,
  2128849065,
  2128925811,
  2128994934,
  2129056501,
  2129110560,
  2129157136,
  2129196237,
  2129227847,
  2129251929,
  2129268426,
  2129277255,
  2129278312,
  2129271467,
  2129256561,
  2129233410,
  2129201800,
  2129161480,
  2129112170,
  2129053545,
  2128985244,
  2128906855,
  2128817916,
  2128717911,
  2128606255,
  2128482298,
  2128345305,
  2128194452,
  2128028813,
  2127847342,
  2127648860,
  2127432031,
  2127195339,
  2126937058,
  2126655214,
  2126347546,
  2126011445,
  2125643893,
  2125241376,
  2124799783,
  2124314271,
  2123779094,
  2123187386,
  2122530867,
  2121799464,
  2120980787,
  2120059418,
  2119015917,
  2117825402,
  2116455471,
  2114863093,
  2112989789,
  2110753906,
  2108037662,
  2104664315,
  2100355223,
  2094642347,
  2086670106,
  2074676188,
  2054300022,
  2010539237,
  1991057938
};

static const double wtab[128] = {
  1.26809284419130464388715e-10,
  1.68975177699894785991205e-10,
  1.98626884400520997542308e-10,
  2.22324317905021774029349e-10,
  2.42449361237299124520446e-10,
  2.60161318991108784172029e-10,
  2.76119887103316364403791e-10,
  2.90739628164604409183491e-10,
  3.04299704132155695474083e-10,
  3.16997952128715914041369e-10,
  3.28980205260965989233567e-10,
  3.40357381208744496806915e-10,
  3.5121602212754547198276e-10,
  3.61625099496985034075789e-10,
  3.7164057634131740989886e-10,
  3.8130856430312512596469e-10,
  3.90667568091865375507411e-10,
  3.99750118692429207049947e-10,
  4.08583986152762403441697e-10,
  4.17193096394761711796244e-10,
  4.2559823533929967155513e-10,
  4.33817597386126318274218e-10,
  4.41867218119051285930235e-10,
  4.49761319620595590627586e-10,
  4.57512588939983576913862e-10,
  4.65132404808254685733004e-10,
  4.72631023842515177845392e-10,
  4.80017734717789846812131e-10,
  4.87300986774535972100332e-10,
  4.94488498048679623648895e-10,
  5.01587346606858846404499e-10,
  5.08604048237462416940879e-10,
  5.15544622914649265825864e-10,
  5.22414651965840750416997e-10,
  5.2921932749593422526623e-10,
  5.35963495326682955529176e-10,
  5.42651692477542398619773e-10,
  5.49288180030165567763533e-10,
  5.558769720717203800613e-10,
  5.62421861294078455665e-10,
  5.68926441730448309229428e-10,
  5.75394129033424507317127e-10,
  5.81828178635022485605483e-10,
  5.88231702077215644979885e-10,
  5.946076817585619713717e-10,
  6.00958984306954383303572e-10,
  6.07288372758972243869223e-10,
  6.13598517701655156445591e-10,
  6.19892007511889795627399e-10,
  6.26171357811294981823082e-10,
  6.32439020239944976690067e-10,
  6.38697390640029809331121e-10,
  6.44948816730244385497167e-10,
  6.51195605343024410182548e-10,
  6.5744002928946178913137e-10,
  6.63684333910635473870468e-10,
  6.69930743369022996079553e-10,
  6.76181466729481002582792e-10,
  6.82438703875893087890625e-10,
  6.88704651306894540121957e-10,
  6.94981507852028779712267e-10,
  7.01271480348217422256114e-10,
  7.0757678931549693487351e-10,
  7.13899674670564047387398e-10,
  7.20242401516765104999874e-10,
  7.26607266049757919725128e-10,
  7.32996601619175505782248e-10,
  7.39412784988247184946175e-10,
  7.45858242835512877750422e-10,
  7.52335458545541839863296e-10,
  7.5884697933899164024717e-10,
  7.65395423796485480355862e-10,
  7.71983489835731473830826e-10,
  7.78613963207161150498644e-10,
  7.85289726580253953964804e-10,
  7.92013769300794668343784e-10,
  7.98789197908768685610908e-10,
  8.05619247517661906297502e-10,
  8.12507294168871132875681e-10,
  8.19456868290077836679449e-10,
  8.26471669404194339043048e-10,
  8.33555582256344606265481e-10,
  8.40712694550887095719094e-10,
  8.47947316519452759044146e-10,
  8.55264002575252225961458e-10,
  8.62667575349606136671056e-10,
  8.70163152455138977864435e-10,
  8.7775617637805139553244e-10,
  8.85452447971476993517132e-10,
  8.93258164105812101766512e-10,
  9.01179960133461375303735e-10,
  9.09224957948964510968321e-10,
  9.17400820576452288875581e-10,
  9.25715814401889557171101e-10,
  9.34178880396749284264771e-10,
  9.42799715964558477625693e-10,
  9.51588869397840228657448e-10,
  9.60557849381101949827052e-10,
  9.69719252543395685084955e-10,
  9.7908691278891585944201e-10,
  9.88676077066822889981155e-10,
  9.98503613451617468033378e-10,
  1.00858825898954679151739e-09,
  1.01895091686026221572322e-09,
  1.029615015198815525789e-09,
  1.04060694369816080679612e-09,
  1.05195658927100222793561e-09,
  1.06369799919131063185578e-09,
  1.0758702101628306115789e-09,
  1.08851829605900251093084e-09,
  1.10169470781180427572374e-09,
  1.11546100955804220517514e-09,
  1.12989016134767401635861e-09,
  1.14506957000510305832519e-09,
  1.16110524260064156767471e-09,
  1.17812756094404804380418e-09,
  1.1962995053835395410798e-09,
  1.2158286983280504138695e-09,
  1.23698562907902174626409e-09,
  1.26013233005941022983024e-09,
  1.28576968441910755421728e-09,
  1.31462018496634738420194e-09,
  1.34778395621975447096936e-09,
  1.38706353150541752603784e-09,
  1.43574031918040183413463e-09,
  1.50086590302112444517283e-09,
  1.60309479380801884523794e-09,
  1.72904052154177951669258e-09
};


// funzione gamma per valori n/2




//genera l'istogramma (intervallo [a,b], numero di intervalli dell'istogramma z, numero di gruppi nc di N numeri da ricevere da *sorgente, array dell'istogramma, funzione generatrice)
void istg(double a, double b, unsigned int z, unsigned int nc,unsigned int N, unsigned int *ist, void (*sorgente)(double *),double *mu, double *var){
unsigned int j,jj;
double k,*x;
*mu=0;
*var=0;
x= (double *) malloc(sizeof(double)*N);
 for (j=0;j<nc;j++){
  sorgente(x);
   for (jj=0;jj<N;jj++){
    *mu=*mu+x[jj];
    *var=*var+x[jj]*x[jj];
    k=((x[jj]-a)/(b-a))*z;
     if (k<0) {k=0;} else if (k>z-1) {k=z-1;}
    ist[(unsigned int) k]++;
   }
 }
*mu=*mu/(double) (nc*N);
*var=*var/(double) (nc*N-1);
*var=(*var-(*mu)*(*mu))*nc*N/(nc*N-1);
free(x);
}

void istg_(double a,/// estremo inferiore dell'intervallo da istogrammare
           double b,///estremo superiore dell'intervallo da istogrammare
           unsigned int z,/// numero di intervalli dell'istogramma
           unsigned int nc,/// numero di dati da istogrammare (numero di chiamate alla funzione sorgente)
           unsigned int *ist,/// array che conterrà dopo l'istogramma calcolato. E' responsabilità dell'utente la sua inizializzazione (a zero o a quello che si vuole)
           double (*sorgente)(), ///funzione che viene chiamata nc volte per leggere i dati da istogrammare
           double *mu,/// qui poi scrivo il valore medio dei dati istogrammati
           double *var/// qui metto la varianza dei dati istogrammati
           ){
unsigned int j,jj;
double k,x;
*mu=0;
*var=0;
 for (j=0;j<nc;j++){
  x=sorgente();
  *mu=*mu+x;
  *var=*var+x*x;
    k=(x-a)/(b-a)*z;
     if (k<0) {k=0;} else if (k>z-1) {k=z-1;}
    ist[(unsigned int) k]++;
 }
*mu=*mu/(double) nc;
*var=*var/(double) (nc-1);
*var=(*var-(*mu)*(*mu))*nc/(nc-1);
}

void istg__(double a, double b, unsigned int z, unsigned int nc, unsigned int *ist, float (*sorgente)(),float *mu, float *var){
unsigned int j,jj;
float k,x;
*mu=0;
*var=0;
 for (j=0;j<nc;j++){
  x=sorgente();
  *mu=*mu+x;
  *var=*var+x*x;
    k=(x-a)/(b-a)*z;
     if (k<0) {k=0;} else if (k>z-1) {k=z-1;}
    ist[(unsigned int) k]++;
 }
*mu=*mu/(float) nc;
*var=*var/(float) (nc-1);
*var=(*var-(*mu)*(*mu))*nc/(nc-1);
}


void istogramma(float a /// estremo inferiore dell'intervallo da istogrammare
                , float b ///estremo superiore dell'intervallo da istogrammare
                , float *xx /// scrive qui gli estremi di ciascun intervallino
                , unsigned int z /// numero di intervalli dell'istogramma
                , unsigned int nc /// numero di dati da istogrammare (lunghezza dell'array sorgente)
                , unsigned int *ist /// array che conterrà dopo l'istogramma calcolato. deve essere inizializzato a zero
                , float *sorgente /// array lungo zc con i dati da istogrammare
                ,float *mu /// qui poi scrivo il valore medio dei dati istogrammati
                , float *var /// qui metto la varianza dei dati istogrammati
                ){
unsigned int j,jj,c=0;
float k,x,delta,mean=0,M2=0,v=0;
*mu=0;
*var=0;
 for (j=0;j<nc;j++){
  x=sorgente[j];
    k=((x-a)/(b-a))*(float)z;
     if (k<0 || k>=z) {fprintf (stderr,"#Valore %f (indice %i) fuori dall'intervallo [%f %f] dell'istogramma (indice %f %i)\n",x,j,a,b,k,z);} else {
      c++;
      ist[(unsigned int) k]++;
        delta = x - mean;
        mean = mean + delta/c;
        M2 = M2 + delta*(x - mean);
     }
 }
 for (j=0;j<z+1;j++){
  xx[j]=a+(b-a)/(double)z*j;
 }
*var=M2/((float)c - 1);
*mu=mean;
}

void covarianza(unsigned int n,float *x, float *y, float *mu, float *var, float * covar){
    unsigned int i;
    mu[0]=0;
    mu[1]=0;
    var[0]=0;
    var[1]=0;
    covar[0]=0;
    for (i=0;i<n;i++){
        mu[0]+=x[i];
        mu[1]+=y[i];
    }
    mu[0]=mu[0]/(float) n;
    mu[1]=mu[1]/(float) n;
    for (i=0;i<n;i++){
        var[0]+=(x[i]-mu[0])*(x[i]-mu[0]);
        var[1]+=(y[i]-mu[1])*(y[i]-mu[1]);
        *covar+=(x[i]-mu[0])*(y[i]-mu[1]);
    }
    var[0]=var[0]/(n-1);
    var[1]=var[1]/(n-1);
    covar[0]=covar[0]/(n-1);
}

// rudimentale generatore, che genera numeri con qualita' non molto alta
uint32_t lcg() {
static uint32_t x_=54321; // static dice di conservare il valore della variabile fra una chiamata e l'altra della funzione
x_=1664525*x_+1013904223;
return x_;
}


// questa funzione inizializza il generatore, quello buono
void init_cmwc4096(){
int i;
    for (i=0;i<4096;i++){
        Q_[i]=SHR3;
    }
}

// questa funzione genera un numero intero pseudo-casuale a 32 bit (compreso fra 0 e 2^32-1) -- algoritmo di Marsaglia
uint32_t cmwc4096(){
uint64_t t, a=18782LL;
static uint32_t i=4095,c_=362436;
uint32_t x,r=0xfffffffe;
    i=(i+1)&4095;
    t=a*Q_[i]+c_;
    c_=(t>>32);
    x=t+c_;
    if(x<c_){
        x++;c_++;
    }
    return(Q_[i]=r-x);
}

// questa funzione combina due numeri interi a 32 bit in un numero intero a 64 bit, usando le operazioni fra bit.
uint64_t rnd64(void){
uint64_t t,r=0;
r=( ( (uint64_t) cmwc4096() << 32 )& 0xFFFFFFFF00000000ull) | ( (uint64_t) cmwc4096() & 0x00000000FFFFFFFFull); // << sposta i bit a sinistra della quantita' intera specificata. | esegue un "or" bit per bit, & "and"
return r;
}

// prende direttamente i bit da un numero intero casuale e li mette dentro la mantissa del numero float, impostando anche l'esponente. Il numero generato e' compreso fra 0.5 e 1 (il numero e' scritto in base 2).
double rnd_double() {
double t;
uint64_t r;
uint32_t rr;
    if (sizeof(double)==8) { // float a 64 bit
        r=rnd64();
//			      001111111110=0x3FE
        r= (r>>12)& 0x000FFFFFFFFFFFFFull  | 0x3FE0000000000000ull;
        t= * (double*) &r;
    }else if (sizeof(double)==4) { // float a 32 bit
        rr=cmwc4096();
//		              001111110000 = 0x3F000000
        rr=(rr>>9)& 0b00000000011111111111111111111111 | 0x3F000000;
        t=*(double*)&rr;
    } else {
        fprintf(stderr,"Errore: tipo double non riconosciuto.\n");
    }
return t;
}

// stessa cosa con il tipo a 32 bit
float rnd_single() {
uint32_t r;
float t;
    if (sizeof(float)==4){
        r=cmwc4096();
//		            001111110000
        r=(r>>9)& 0b00000000011111111111111111111111 | 0x3F000000;
        t=*(float*)&r;
    } else {
        fprintf(stderr,"Errore: tipo float non riconosciuto.\n");
    }
return t;
}

// prima rudimentale implementazione della trasformazione di box-muller
void rnd_gauss(double *ris) {
double u,v,r,t;
    while(1) { // scegliamo due numeri non contemporaneamente nulli u,v in [-1,1] tali che (u^2+v^2<1)
        u=(rnd_double()-0.75)*4; // come sono distribuiti i numeri casuali nell'intervallo [-1,1]?
        v=(rnd_double()-0.75)*4;
        r= u*u+v*v;
            if (r>0 && r<1) break;
    }
    t=sqrt(-2.0*log(r)/r); // formula della trasformazione.
    ris[0]=u*t;
    ris[1]=v*t;
}

/* implementazione avanzata di un generatore di numeri pseudo-casuali con distribuzione gaussiana (http://www.jstatsoft.org/v05/i08/paper (George Marsaglia; Wai Wan Tsang (2000). "The Ziggurat Method for Generating Random Variables". Journal of Statistical Software 5 (8).)
*/

double normal_gauss(){
int i,indice;
double x,y;

i=cmwc4096();
indice=i&127;
x=wtab[indice]*i;
    if (abs(i)<ktab[indice]){ // sono completamente dentro la gaussiana, accetto subito!
        return i*wtab[indice];
    }else if (i==127){ // aiai sono finito sulla coda
            do {
                x=-log(UNI)/PARAM_R;
                y=-log(UNI);
            } while (y+y<x*x);
            return (i>0)? PARAM_R+x : -PARAM_R-x;
    } else if (ytab[indice]+UNI*(ytab[indice+1]-ytab[indice]) < exp(-.5*x*x) ) { // sono sul bordo esterno del rettangolino, devo vedere dove cade la y per decidere se sono dentro la gaussiana o no!
        return x;
    }else {
        return normal_gauss();
    }

}

void ols2 (unsigned int n, float *x, float *y, float *ey, float *m, float *q, float *em, float *eq, float *chi2){
unsigned int j;
float s01=0,s10=0,s00=0,s11=0,s20=0;
for (j=0;j<n;j++){
    s01+=y[j]/ey[j];
    s10+=x[j]/ey[j];
    s00+=1/ey[j];
    s11+=x[j]*y[j]/ey[j];
    s20+=x[j]*x[j]/ey[j];
}
 *m=(s11*s00-s10*s01)/(s20*s00-s10*s10);
 *em=s00/(s20*s00-s10*s10);
 *q=(s20*s01-s11*s10)/(s20*s00-s10*s10);
 *eq=s20/(s20*s00-s10*s10);
 *chi2=0;
for (j=0;j<n;j++){
    *chi2+=(y[j]-*m*x[j]-*q)*(y[j]-*m*x[j]-*q)/(ey[j]*ey[j]);
}
}

void ols (unsigned int n, double *x, double *y, double *m, double *q, double *em, double *eq, double *chi2){
unsigned int j;
double s01=0,s10=0,s00=0,s11=0,s20=0;
for (j=0;j<n;j++){
    s01+=y[j];
    s10+=x[j];
    s11+=x[j]*y[j];
    s20+=x[j]*x[j];
}
 s00=n;
 *m=(s11*s00-s10*s01)/(s20*s00-s10*s10);
 *em=s00/(s20*s00-s10*s10);
 *q=(s20*s01-s11*s10)/(s20*s00-s10*s10);
 *eq=s20/(s20*s00-s10*s10);
 *chi2=0;
for (j=0;j<n;j++){
    *chi2+=(y[j]-*m*x[j]-*q)*(y[j]-*m*x[j]-*q);
}
}

// generico integrale "pesato"

float intg(float(*f)(float),float a, float b, unsigned int n, float *pesi, unsigned int p) {
unsigned int i,ii,nn;
float h,integ=0;
nn=n/(p-1);
h=(b-a)/(p-1)/nn;
    for (i=0;i<nn;i++) {
        for (ii=0;ii<p;ii++) {
            integ+=f(a+h*(i*(p-1)+ii))*pesi[ii];
        }
    }
integ*=h;
return integ;
}



/**Generico integrale pesato (con le formule di Newton-Cotes)
 *
 * Implementa anche un generico puntatore che può essere usato per passare informazioni alla funzione integranda.
 * In questo programma questo serve a passare il puntatore della classe della funzione da integrare e quindi poter
 * integrare anche una funzione membro di una classe C++ (sporco trucco).
 * */
float intg_C(float(*f)(float,void*), ///< funzione da integrare
             float a, ///< estremo più basso
             float b, ///< estremo più alto
             unsigned int n, ///< numero di intervallini da usare per fare l'integrale
             float *pesi, ///< array con i pesi da usare per l'integrale
             unsigned int p, ///< numero di pesi
             void * C ///< puntatore da passare alla funzione da integrare
             ) {
unsigned int i,ii,nn;
float h,integ=0;
nn=n/(p-1);
h=(b-a)/(p-1)/nn;
    for (i=0;i<nn;i++) {
        for (ii=0;ii<p;ii++) {
            integ+=f(a+h*(i*(p-1)+ii),C)*pesi[ii];
        }
    }
integ*=h;
return integ;
}

#ifndef ANALISI
/**Calcola il chi^2 dato un istogramma e un distribuzione da testare*/
void pearson(unsigned int n, ///< numero di bin dell'istogramma
             float *x, ///< estremi di ogni bin (n+1 elementi)
             unsigned int *ist, ///< array dell'istogramma
             float *chi2, ///< qui memorizza il chi^2
             unsigned int *nu, ///< qui memorizza i gradi di libertà del chi^2 (se i bin hanno pochi elementi l'algoritmo li conta assieme ad altri)
             float (*fdistrib)(float) ///< puntatore alla distribuzione da testare
             ){
unsigned int i=0,cont=0,cont_inizio=0,cont_inizioP,contP,ne=0;
float p,chi2P=0,chi2Pold,pt=0;
*chi2=0;
*nu=0;
//printf ("#passo	x		p		ev.mis.	ev.calc		chi2+		ampiezza");
    for (i=0;i<n;i++) ne+=ist[i];
    for (i=0;i<n;i++){
        cont+=ist[i];
        if (cont>=10){ // almeno 10 eventi per intervallo
            *nu=*nu+1;
            p=intg(fdistrib,x[cont_inizio],x[i+1],50,nc4,4);
            pt+=p;
//			printf("\n%i	%1.6f	%1.6f	%i	%1.6f	%1.6f	%1.6f",i,x[cont_inizio],p,cont,p*ne,(cont-p*ne)*(cont-p*ne)/(p*ne),x[i+1]-x[cont_inizio]);
            chi2Pold=chi2P; // nel caso negli ultimi intervalli ci siano meno di 10 eventi fondo gli ultimi intervalli, e il chi2 di base diventa questo
            chi2P+=(cont-p*ne)*(cont-p*ne)/(p*ne);
            contP=cont; // nel caso negli ultimi intervalli ci siano meno di 10
            cont=0;
            cont_inizioP=cont_inizio; // nel caso ultimi meno di 10
            cont_inizio=i+1;
        }
    }
    if (cont!=0 && cont < 10){ // negli ultimi intervalli ci sono meno di 10 elementi
        p=intg(fdistrib,x[cont_inizioP],x[n],50,nc4,4);
//		printf("\n!la riga precedente e' contata nella seguente!\n%i	%1.6f	%1.6f	%i	%1.6f	%1.6f	%1.6f\n",i,x[cont_inizioP],p,cont+contP,p*ne,(cont+contP-p*ne)*(cont+contP-p*ne)/(p*ne),x[n]-x[cont_inizioP]);
        *chi2=chi2Pold+(cont+contP-p*ne)*(cont+contP-p*ne)/(p*ne);
    } else {
//		printf("\n");
        *chi2=chi2P;
    }
//printf("#Somma probabilita' %f,chi2 = %f, chi2P = %f\n",pt,*chi2,chi2P);
}

// c++ class function pointer
void pearson_C(unsigned int n, float *x, unsigned int *ist,float *chi2,unsigned int *nu, float (*fdistrib)(float,void*),void* C){
unsigned int i=0,cont=0,cont_inizio=0,cont_inizioP,contP,ne=0;
float p,chi2P=0,chi2Pold,pt=0;
*chi2=0;
*nu=0;
//printf ("#passo	x		p		ev.mis.	ev.calc		chi2+		ampiezza");
    for (i=0;i<n;i++) ne+=ist[i];
    for (i=0;i<n;i++){
        cont+=ist[i];
        if (cont>=10){ // almeno 10 eventi per intervallo
            *nu=*nu+1;
            p=intg_C(fdistrib,x[cont_inizio],x[i+1],50,nc4,4,C);
            pt+=p;
//			printf("\n%i	%1.6f	%1.6f	%i	%1.6f	%1.6f	%1.6f",i,x[cont_inizio],p,cont,p*ne,(cont-p*ne)*(cont-p*ne)/(p*ne),x[i+1]-x[cont_inizio]);
            chi2Pold=chi2P; // nel caso negli ultimi intervalli ci siano meno di 10 eventi fondo gli ultimi intervalli, e il chi2 di base diventa questo
            chi2P+=(cont-p*ne)*(cont-p*ne)/(p*ne);
            contP=cont; // nel caso negli ultimi intervalli ci siano meno di 10
            cont=0;
            cont_inizioP=cont_inizio; // nel caso ultimi meno di 10
            cont_inizio=i+1;
        }
    }
    if (cont!=0 && cont < 10){ // negli ultimi intervalli ci sono meno di 10 elementi
        p=intg_C(fdistrib,x[cont_inizioP],x[n],50,nc4,4,C);
//		printf("\n!la riga precedente e' contata nella seguente!\n%i	%1.6f	%1.6f	%i	%1.6f	%1.6f	%1.6f\n",i,x[cont_inizioP],p,cont+contP,p*ne,(cont+contP-p*ne)*(cont+contP-p*ne)/(p*ne),x[n]-x[cont_inizioP]);
        *chi2=chi2Pold+(cont+contP-p*ne)*(cont+contP-p*ne)/(p*ne);
    } else {
//		printf("\n");
        *chi2=chi2P;
    }
//printf("#Somma probabilita' %f,chi2 = %f, chi2P = %f\n",pt,*chi2,chi2P);
}


/* n dispari
   __              ____
  |   / n \       /___ |  (n - 2)!!
  |   |---|  =  \/ | |   -----------
  |   \ 2 /                (n-1)/2
                          2

*/


void gamma_n2 (unsigned int n, mpfr_t gamma) {
mpz_t t1;
mpfr_t f1,f2;

mpz_init(t1);
mpfr_inits(f1,f2,NULL);
    if (n%2==0) { // se n e' pari restituisce il fattoriale
        tr+=mpfr_fac_ui(gamma,n/2-1,MPFR_RNDN);
    } else { // se n e' dispari usa la formula molto più lento. Sicuramente si può fare meglio!
        // 2^{(n-1)/2}
        tr+=mpfr_set_ui(f2,(n-1),MPFR_RNDN);
        tr+=mpfr_div_ui(f2,f2,2,MPFR_RNDN);
        tr+=mpfr_exp2(f1,f2,MPFR_RNDN);
        // (n-2)!!
        mpz_2fac_ui(t1,n-2);
        tr+=mpfr_set_z(f2,t1,MPFR_RNDN);
        // f2/f1
        tr+=mpfr_div(f1,f2,f1,MPFR_RNDN);
        // sqrt(pi)
        tr+=mpfr_set_si(f2,-1,MPFR_RNDN);
        tr+=mpfr_acos(f2,f2,MPFR_RNDN);
        tr+=mpfr_sqrt(f2,f2,MPFR_RNDN);
        // f1*f2
        tr+=mpfr_mul(gamma,f1,f2,MPFR_RNDN);
    }
}

/**Calcola i quantili della distribuzione chi quadro
 *
 * Usa un metodo brutale: va avanti a integrare sempre di più fino a quando il valore integrato non supera
 * il valore del quantile richiesto. Quindi restituisce il valore dell'ascissa a cui era arrivato.
 * */

float chi2_quantile(float p, ///< quantile da calcolare
                    unsigned int nu ///< gradi di libertà della distribuzione chi quadro
                    ) {

unsigned int i,ii,nn,n=200000;

if (nu<201 && p>0.9499 && p<0.9501) return chi2_095q[nu];
if (nu<201 && p>0.5999 && p<0.6001) return chi2_060q[nu];
if (nu<201 && p>0.4999 && p<0.5001) return chi2_050q[nu];
if (nu<201 && p>0.3999 && p<0.4001) return chi2_040q[nu];

mpfr_t c,t1,t2,t3,h,x,integ,b,a,n4[4];
mpfr_set_default_prec(PRECISIONE);
mpfr_inits(c,t1,t2,t3,h,x,integ,b,a,NULL);
for(i=0;i<4;i++) mpfr_init(n4[i]);

tr+=mpfr_set_ui(b,100,MPFR_RNDN);
tr+=mpfr_set_ui(a,0,MPFR_RNDN);
tr+=mpfr_set_ui(integ,0,MPFR_RNDN);

//c=(exp2((float)nu/2)*gamma_n2(nu));
tr+=mpfr_set_ui(t1,nu,MPFR_RNDN);
tr+=mpfr_div_ui(t1,t1,2,MPFR_RNDN);
tr+=mpfr_exp2(t1,t1,MPFR_RNDN);
gamma_n2(nu,t2);
tr+=mpfr_mul(c,t1,t2,MPFR_RNDN);
//static float nc4[4] = {3.0/8.0, 9.0/8.0, 9.0/8.0, 3.0/8.0};
tr+=mpfr_set_ui(n4[0],3,MPFR_RNDN);
tr+=mpfr_set_ui(n4[1],9,MPFR_RNDN);
tr+=mpfr_set_ui(n4[2],9,MPFR_RNDN);
tr+=mpfr_set_ui(n4[3],3,MPFR_RNDN);
for(i=0;i<4;i++) tr+=mpfr_div_ui(n4[i],n4[i],8,MPFR_RNDN);

while (1) {
    nn=n/3;
//	h=(b-a)/(3)/nn;
//	integ=0;
    tr+=mpfr_sub(h,b,a,MPFR_RNDN);
    tr+=mpfr_div_ui(h,h,3*nn,MPFR_RNDN);
    for (i=0;i<nn;i++) {
        tr+=mpfr_set_ui(t3,0,MPFR_RNDN);
        for (ii=0;ii<4;ii++) {
//			x = a+h*(i*(3)+ii);
            tr+=mpfr_mul_ui(x,h,i*3+ii,MPFR_RNDN);
            tr+=mpfr_add(x,x,a,MPFR_RNDN);
//			integ+=exp(-x/2)*pow(x,(float)nu/2-1)*nc4[ii]*h/c;
            tr+=mpfr_div_ui(t1,x,2,MPFR_RNDN);
            tr+=mpfr_neg(t1,t1,MPFR_RNDN);
            tr+=mpfr_exp(t1,t1,MPFR_RNDN);
            tr+=mpfr_set_ui(t2,nu,MPFR_RNDN);
            tr+=mpfr_div_ui(t2,t2,2,MPFR_RNDN);
            tr+=mpfr_sub_ui(t2,t2,1,MPFR_RNDN);
            tr+=mpfr_pow(t2,x,t2,MPFR_RNDN);
            tr+=mpfr_mul(t1,t1,t2,MPFR_RNDN);
            tr+=mpfr_mul(t1,t1,n4[ii],MPFR_RNDN);
            tr+=mpfr_add(t3,t3,t1,MPFR_RNDN);
        }
            tr+=mpfr_mul(t3,t3,h,MPFR_RNDN);
            tr+=mpfr_div(t3,t3,c,MPFR_RNDN);
            tr+=mpfr_add(integ,integ,t3,MPFR_RNDN);
        if (mpfr_cmp_d(integ,(double)p)>=0) {return mpfr_get_flt(x,MPFR_RNDN);}
    }
//	b*=2;
    tr+=mpfr_add_ui(b,b,100,MPFR_RNDN);
    tr+=mpfr_add_ui(a,a,100,MPFR_RNDN); // se non ho ancora trovato il valore procedo
//	n*=2;
}
}
//ifndef ANALISI
#endif
