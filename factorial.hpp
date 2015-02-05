#ifndef FACTORIAL_HPP
#define FACTORIAL_HPP

#include <cassert>

// Table adopted and modified from GNU Scientific Library Version 1.9 (gamma.c).

const int FACTORIAL_NMAX(170);
static const double factorials[FACTORIAL_NMAX + 1] =
{
    1.0,
    1.0,
    2.0,
    6.0,
    24.0,
    120.0,
    720.0,
    5040.0,
    40320.0,
    362880.0,
    3628800.0,
    39916800.0,
    479001600.0,
    6227020800.0,
    87178291200.0,
    1307674368000.0,
    20922789888000.0,
    355687428096000.0,
    6402373705728000.0,
    121645100408832000.0,
    2432902008176640000.0,
    51090942171709440000.0,
    1124000727777607680000.0,
    25852016738884976640000.0,
    620448401733239439360000.0,
    15511210043330985984000000.0,
    403291461126605635584000000.0,
    10888869450418352160768000000.0,
    304888344611713860501504000000.0,
    8841761993739701954543616000000.0,
    265252859812191058636308480000000.0,
    8222838654177922817725562880000000.0,
    263130836933693530167218012160000000.0,
    8683317618811886495518194401280000000.0,
    2.95232799039604140847618609644e38,
    1.03331479663861449296666513375e40,
    3.71993326789901217467999448151e41,
    1.37637530912263450463159795816e43,
    5.23022617466601111760007224100e44,
    2.03978820811974433586402817399e46,
    8.15915283247897734345611269600e47,
    3.34525266131638071081700620534e49,
    1.40500611775287989854314260624e51,
    6.04152630633738356373551320685e52,
    2.65827157478844876804362581101e54,
    1.19622220865480194561963161496e56,
    5.50262215981208894985030542880e57,
    2.58623241511168180642964355154e59,
    1.24139155925360726708622890474e61,
    6.08281864034267560872252163321e62,
    3.04140932017133780436126081661e64,
    1.55111875328738228022424301647e66,
    8.06581751709438785716606368564e67,
    4.27488328406002556429801375339e69,
    2.30843697339241380472092742683e71,
    1.26964033536582759259651008476e73,
    7.10998587804863451854045647464e74,
    4.05269195048772167556806019054e76,
    2.35056133128287857182947491052e78,
    1.38683118545689835737939019720e80,
    8.32098711274139014427634118320e81,
    5.07580213877224798800856812177e83,
    3.14699732603879375256531223550e85,
    1.982608315404440064116146708360e87,
    1.268869321858841641034333893350e89,
    8.247650592082470666723170306800e90,
    5.443449390774430640037292402480e92,
    3.647111091818868528824985909660e94,
    2.480035542436830599600990418570e96,
    1.711224524281413113724683388810e98,
    1.197857166996989179607278372170e100,
    8.504785885678623175211676442400e101,
    6.123445837688608686152407038530e103,
    4.470115461512684340891257138130e105,
    3.307885441519386412259530282210e107,
    2.480914081139539809194647711660e109,
    1.885494701666050254987932260860e111,
    1.451830920282858696340707840860e113,
    1.132428117820629783145752115870e115,
    8.946182130782975286851441715400e116,
    7.156945704626380229481153372320e118,
    5.797126020747367985879734231580e120,
    4.753643337012841748421382069890e122,
    3.945523969720658651189747118010e124,
    3.314240134565353266999387579130e126,
    2.817104114380550276949479442260e128,
    2.422709538367273238176552320340e130,
    2.107757298379527717213600518700e132,
    1.854826422573984391147968456460e134,
    1.650795516090846108121691926250e136,
    1.485715964481761497309522733620e138,
    1.352001527678402962551665687590e140,
    1.243841405464130725547532432590e142,
    1.156772507081641574759205162310e144,
    1.087366156656743080273652852570e146,
    1.032997848823905926259970209940e148,
    9.916779348709496892095714015400e149,
    9.619275968248211985332842594960e151,
    9.426890448883247745626185743100e153,
    9.332621544394415268169923885600e155,
    9.33262154439441526816992388563e157,
    9.42594775983835942085162312450e159,
    9.61446671503512660926865558700e161,
    9.90290071648618040754671525458e163,
    1.02990167451456276238485838648e166,
    1.08139675824029090050410130580e168,
    1.146280563734708354534347384148e170,
    1.226520203196137939351751701040e172,
    1.324641819451828974499891837120e174,
    1.443859583202493582204882102460e176,
    1.588245541522742940425370312710e178,
    1.762952551090244663872161047110e180,
    1.974506857221074023536820372760e182,
    2.231192748659813646596607021220e184,
    2.543559733472187557120132004190e186,
    2.925093693493015690688151804820e188,
    3.393108684451898201198256093590e190,
    3.96993716080872089540195962950e192,
    4.68452584975429065657431236281e194,
    5.57458576120760588132343171174e196,
    6.68950291344912705758811805409e198,
    8.09429852527344373968162284545e200,
    9.87504420083360136241157987140e202,
    1.21463043670253296757662432419e205,
    1.50614174151114087979501416199e207,
    1.88267717688892609974376770249e209,
    2.37217324288004688567714730514e211,
    3.01266001845765954480997707753e213,
    3.85620482362580421735677065923e215,
    4.97450422247728744039023415041e217,
    6.46685548922047367250730439554e219,
    8.47158069087882051098456875820e221,
    1.11824865119600430744996307608e224,
    1.48727070609068572890845089118e226,
    1.99294274616151887673732419418e228,
    2.69047270731805048359538766215e230,
    3.65904288195254865768972722052e232,
    5.01288874827499166103492629211e234,
    6.91778647261948849222819828311e236,
    9.61572319694108900419719561353e238,
    1.34620124757175246058760738589e241,
    1.89814375907617096942852641411e243,
    2.69536413788816277658850750804e245,
    3.85437071718007277052156573649e247,
    5.55029383273930478955105466055e249,
    8.04792605747199194484902925780e251,
    1.17499720439091082394795827164e254,
    1.72724589045463891120349865931e256,
    2.55632391787286558858117801578e258,
    3.80892263763056972698595524351e260,
    5.71338395644585459047893286526e262,
    8.62720977423324043162318862650e264,
    1.31133588568345254560672467123e267,
    2.00634390509568239477828874699e269,
    3.08976961384735088795856467036e271,
    4.78914290146339387633577523906e273,
    7.47106292628289444708380937294e275,
    1.17295687942641442819215807155e278,
    1.85327186949373479654360975305e280,
    2.94670227249503832650433950735e282,
    4.71472363599206132240694321176e284,
    7.59070505394721872907517857094e286,
    1.22969421873944943411017892849e289,
    2.00440157654530257759959165344e291,
    3.28721858553429622726333031164e293,
    5.42391066613158877498449501421e295,
    9.00369170577843736647426172359e297,
    1.50361651486499904020120170784e300,
    2.52607574497319838753801886917e302,
    4.26906800900470527493925188890e304,
    7.25741561530799896739672821113e306,
    /*
      1.24101807021766782342484052410e309,
      2.13455108077438865629072570146e311,
      3.69277336973969237538295546352e313,
      6.42542566334706473316634250653e315,
      1.12444949108573632830410993864e318,
      1.97903110431089593781523349201e320,
      3.50288505463028580993296328086e322,
      6.23513539724190874168067463993e324,
      1.11608923610630166476084076055e327,
      2.00896062499134299656951336898e329,
      3.63621873123433082379081919786e331,
      6.61791809084648209929929094011e333,
      1.21107901062490622417177024204e336,
      2.22838537954982745247605724535e338,
      4.12251295216718078708070590390e340,
      7.66787409103095626397011298130e342,
      1.43389245502278882136241112750e345,
      2.69571781544284298416133291969e347,
      5.09490667118697324006491921822e349,
      9.68032267525524915612334651460e351,
      1.84894163097375258881955918429e354,
      3.54996793146960497053355363384e356,
      6.85143810773633759312975851330e358,
      1.32917899290084949306717315158e361,
      2.59189903615665651148098764559e363,
      5.08012211086704676250273578535e365,
      1.00078405584080821221303894971e368,
      1.98155243056480026018181712043e370,
      3.94328933682395251776181606966e372,
      7.88657867364790503552363213932e374,
      */
};

static const double factorials_r[FACTORIAL_NMAX + 1] =
{
    1.0 / 1.0,
    1.0 / 1.0,
    1.0 / 2.0,
    1.0 / 6.0,
    1.0 / 24.0,
    1.0 / 120.0,
    1.0 / 720.0,
    1.0 / 5040.0,
    1.0 / 40320.0,
    1.0 / 362880.0,
    1.0 / 3628800.0,
    1.0 / 39916800.0,
    1.0 / 479001600.0,
    1.0 / 6227020800.0,
    1.0 / 87178291200.0,
    1.0 / 1307674368000.0,
    1.0 / 20922789888000.0,
    1.0 / 355687428096000.0,
    1.0 / 6402373705728000.0,
    1.0 / 121645100408832000.0,
    1.0 / 2432902008176640000.0,
    1.0 / 51090942171709440000.0,
    1.0 / 1124000727777607680000.0,
    1.0 / 25852016738884976640000.0,
    1.0 / 620448401733239439360000.0,
    1.0 / 15511210043330985984000000.0,
    1.0 / 403291461126605635584000000.0,
    1.0 / 10888869450418352160768000000.0,
    1.0 / 304888344611713860501504000000.0,
    1.0 / 8841761993739701954543616000000.0,
    1.0 / 265252859812191058636308480000000.0,
    1.0 / 8222838654177922817725562880000000.0,
    1.0 / 263130836933693530167218012160000000.0,
    1.0 / 8683317618811886495518194401280000000.0,
    1.0 / 2.95232799039604140847618609644e38,
    1.0 / 1.03331479663861449296666513375e40,
    1.0 / 3.71993326789901217467999448151e41,
    1.0 / 1.37637530912263450463159795816e43,
    1.0 / 5.23022617466601111760007224100e44,
    1.0 / 2.03978820811974433586402817399e46,
    1.0 / 8.15915283247897734345611269600e47,
    1.0 / 3.34525266131638071081700620534e49,
    1.0 / 1.40500611775287989854314260624e51,
    1.0 / 6.04152630633738356373551320685e52,
    1.0 / 2.65827157478844876804362581101e54,
    1.0 / 1.19622220865480194561963161496e56,
    1.0 / 5.50262215981208894985030542880e57,
    1.0 / 2.58623241511168180642964355154e59,
    1.0 / 1.24139155925360726708622890474e61,
    1.0 / 6.08281864034267560872252163321e62,
    1.0 / 3.04140932017133780436126081661e64,
    1.0 / 1.55111875328738228022424301647e66,
    1.0 / 8.06581751709438785716606368564e67,
    1.0 / 4.27488328406002556429801375339e69,
    1.0 / 2.30843697339241380472092742683e71,
    1.0 / 1.26964033536582759259651008476e73,
    1.0 / 7.10998587804863451854045647464e74,
    1.0 / 4.05269195048772167556806019054e76,
    1.0 / 2.35056133128287857182947491052e78,
    1.0 / 1.38683118545689835737939019720e80,
    1.0 / 8.32098711274139014427634118320e81,
    1.0 / 5.07580213877224798800856812177e83,
    1.0 / 3.14699732603879375256531223550e85,
    1.0 / 1.982608315404440064116146708360e87,
    1.0 / 1.268869321858841641034333893350e89,
    1.0 / 8.247650592082470666723170306800e90,
    1.0 / 5.443449390774430640037292402480e92,
    1.0 / 3.647111091818868528824985909660e94,
    1.0 / 2.480035542436830599600990418570e96,
    1.0 / 1.711224524281413113724683388810e98,
    1.0 / 1.197857166996989179607278372170e100,
    1.0 / 8.504785885678623175211676442400e101,
    1.0 / 6.123445837688608686152407038530e103,
    1.0 / 4.470115461512684340891257138130e105,
    1.0 / 3.307885441519386412259530282210e107,
    1.0 / 2.480914081139539809194647711660e109,
    1.0 / 1.885494701666050254987932260860e111,
    1.0 / 1.451830920282858696340707840860e113,
    1.0 / 1.132428117820629783145752115870e115,
    1.0 / 8.946182130782975286851441715400e116,
    1.0 / 7.156945704626380229481153372320e118,
    1.0 / 5.797126020747367985879734231580e120,
    1.0 / 4.753643337012841748421382069890e122,
    1.0 / 3.945523969720658651189747118010e124,
    1.0 / 3.314240134565353266999387579130e126,
    1.0 / 2.817104114380550276949479442260e128,
    1.0 / 2.422709538367273238176552320340e130,
    1.0 / 2.107757298379527717213600518700e132,
    1.0 / 1.854826422573984391147968456460e134,
    1.0 / 1.650795516090846108121691926250e136,
    1.0 / 1.485715964481761497309522733620e138,
    1.0 / 1.352001527678402962551665687590e140,
    1.0 / 1.243841405464130725547532432590e142,
    1.0 / 1.156772507081641574759205162310e144,
    1.0 / 1.087366156656743080273652852570e146,
    1.0 / 1.032997848823905926259970209940e148,
    1.0 / 9.916779348709496892095714015400e149,
    1.0 / 9.619275968248211985332842594960e151,
    1.0 / 9.426890448883247745626185743100e153,
    1.0 / 9.332621544394415268169923885600e155,
    1.0 / 9.33262154439441526816992388563e157,
    1.0 / 9.42594775983835942085162312450e159,
    1.0 / 9.61446671503512660926865558700e161,
    1.0 / 9.90290071648618040754671525458e163,
    1.0 / 1.02990167451456276238485838648e166,
    1.0 / 1.08139675824029090050410130580e168,
    1.0 / 1.146280563734708354534347384148e170,
    1.0 / 1.226520203196137939351751701040e172,
    1.0 / 1.324641819451828974499891837120e174,
    1.0 / 1.443859583202493582204882102460e176,
    1.0 / 1.588245541522742940425370312710e178,
    1.0 / 1.762952551090244663872161047110e180,
    1.0 / 1.974506857221074023536820372760e182,
    1.0 / 2.231192748659813646596607021220e184,
    1.0 / 2.543559733472187557120132004190e186,
    1.0 / 2.925093693493015690688151804820e188,
    1.0 / 3.393108684451898201198256093590e190,
    1.0 / 3.96993716080872089540195962950e192,
    1.0 / 4.68452584975429065657431236281e194,
    1.0 / 5.57458576120760588132343171174e196,
    1.0 / 6.68950291344912705758811805409e198,
    1.0 / 8.09429852527344373968162284545e200,
    1.0 / 9.87504420083360136241157987140e202,
    1.0 / 1.21463043670253296757662432419e205,
    1.0 / 1.50614174151114087979501416199e207,
    1.0 / 1.88267717688892609974376770249e209,
    1.0 / 2.37217324288004688567714730514e211,
    1.0 / 3.01266001845765954480997707753e213,
    1.0 / 3.85620482362580421735677065923e215,
    1.0 / 4.97450422247728744039023415041e217,
    1.0 / 6.46685548922047367250730439554e219,
    1.0 / 8.47158069087882051098456875820e221,
    1.0 / 1.11824865119600430744996307608e224,
    1.0 / 1.48727070609068572890845089118e226,
    1.0 / 1.99294274616151887673732419418e228,
    1.0 / 2.69047270731805048359538766215e230,
    1.0 / 3.65904288195254865768972722052e232,
    1.0 / 5.01288874827499166103492629211e234,
    1.0 / 6.91778647261948849222819828311e236,
    1.0 / 9.61572319694108900419719561353e238,
    1.0 / 1.34620124757175246058760738589e241,
    1.0 / 1.89814375907617096942852641411e243,
    1.0 / 2.69536413788816277658850750804e245,
    1.0 / 3.85437071718007277052156573649e247,
    1.0 / 5.55029383273930478955105466055e249,
    1.0 / 8.04792605747199194484902925780e251,
    1.0 / 1.17499720439091082394795827164e254,
    1.0 / 1.72724589045463891120349865931e256,
    1.0 / 2.55632391787286558858117801578e258,
    1.0 / 3.80892263763056972698595524351e260,
    1.0 / 5.71338395644585459047893286526e262,
    1.0 / 8.62720977423324043162318862650e264,
    1.0 / 1.31133588568345254560672467123e267,
    1.0 / 2.00634390509568239477828874699e269,
    1.0 / 3.08976961384735088795856467036e271,
    1.0 / 4.78914290146339387633577523906e273,
    1.0 / 7.47106292628289444708380937294e275,
    1.0 / 1.17295687942641442819215807155e278,
    1.0 / 1.85327186949373479654360975305e280,
    1.0 / 2.94670227249503832650433950735e282,
    1.0 / 4.71472363599206132240694321176e284,
    1.0 / 7.59070505394721872907517857094e286,
    1.0 / 1.22969421873944943411017892849e289,
    1.0 / 2.00440157654530257759959165344e291,
    1.0 / 3.28721858553429622726333031164e293,
    1.0 / 5.42391066613158877498449501421e295,
    1.0 / 9.00369170577843736647426172359e297,
    1.0 / 1.50361651486499904020120170784e300,
    1.0 / 2.52607574497319838753801886917e302,
    1.0 / 4.26906800900470527493925188890e304,
    1.0 / 7.25741561530799896739672821113e306,
    /*
      1.0 /   1.24101807021766782342484052410e309,
      1.0 /   2.13455108077438865629072570146e311,
      1.0 /   3.69277336973969237538295546352e313,
      1.0 /   6.42542566334706473316634250653e315,
      1.0 /   1.12444949108573632830410993864e318,
      1.0 /   1.97903110431089593781523349201e320,
      1.0 /   3.50288505463028580993296328086e322,
      1.0 /   6.23513539724190874168067463993e324,
      1.0 /   1.11608923610630166476084076055e327,
      1.0 /   2.00896062499134299656951336898e329,
      1.0 /   3.63621873123433082379081919786e331,
      1.0 /   6.61791809084648209929929094011e333,
      1.0 /   1.21107901062490622417177024204e336,
      1.0 /   2.22838537954982745247605724535e338,
      1.0 /   4.12251295216718078708070590390e340,
      1.0 /   7.66787409103095626397011298130e342,
      1.0 /   1.43389245502278882136241112750e345,
      1.0 /   2.69571781544284298416133291969e347,
      1.0 /   5.09490667118697324006491921822e349,
      1.0 /   9.68032267525524915612334651460e351,
      1.0 /   1.84894163097375258881955918429e354,
      1.0 /   3.54996793146960497053355363384e356,
      1.0 /   6.85143810773633759312975851330e358,
      1.0 /   1.32917899290084949306717315158e361,
      1.0 /   2.59189903615665651148098764559e363,
      1.0 /   5.08012211086704676250273578535e365,
      1.0 /   1.00078405584080821221303894971e368,
      1.0 /   1.98155243056480026018181712043e370,
      1.0 /   3.94328933682395251776181606966e372,
      1.0 /   7.88657867364790503552363213932e374,
      */
};

inline const double factorial(const unsigned int n)
{
    assert(n <= FACTORIAL_NMAX);
    return factorials[n];
}

inline const double factorial_r(const unsigned int n)
{
    assert(n <= FACTORIAL_NMAX);
    return factorials_r[n];
}

#endif /* FACTORIAL_HPP */
