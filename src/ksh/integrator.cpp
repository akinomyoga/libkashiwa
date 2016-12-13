#include "integrator.h"
namespace kashiwa{
//NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
// Supported orders of the Gauss-Legendre quadratures are given in the following table:
//   GaussLegendre<12>
//   GaussLegendre<38>
//   GaussLegendre<100>
//   GaussLegendre<200>
//   GaussLegendre<8>
//   GaussLegendre<16>
//   GaussLegendre<32>
//   GaussLegendre<64>
//   GaussLegendre<128>
//   GaussLegendre<256>

namespace integrator_detail{
  template<>
  point_weight_pair GaussLegendre<12>::data[6]={
    {0.12523340851146892e00, 0.24914704581340279e00 },
    {0.36783149899818019e00, 0.23349253653835481e00 },
    {0.58731795428661745e00, 0.20316742672306592e00 },
    {0.76990267419430469e00, 0.16007832854334623e00 },
    {0.90411725637047486e00, 0.10693932599531843e00 },
    {0.98156063424671925e00, 0.47175336386511827e-01},
  };
  template<>
  point_weight_pair GaussLegendre<38>::data[19]={
    {4.078514790457857E-002,8.152502928038517E-002},
    {1.220840253378676E-001,8.098249377059712E-002},
    {2.025704538921172E-001,7.990103324352750E-002},
    {2.817088097901658E-001,7.828784465821101E-002},
    {3.589724404794348E-001,7.615366354844674E-002},
    {4.338471694323767E-001,7.351269258474258E-002},
    {5.058347179279313E-001,7.038250706689798E-002},
    {5.744560210478074E-001,6.678393797914166E-002},
    {6.392544158296818E-001,6.274093339213346E-002},
    {6.997986803791845E-001,5.828039914699768E-002},
    {7.556859037539705E-001,5.343201991033197E-002},
    {8.065441676053168E-001,4.822806186075808E-002},
    {8.520350219323623E-001,4.270315850467386E-002},
    {8.918557390046323E-001,3.689408159402410E-002},
    {9.257413320485842E-001,3.083950054517547E-002},
    {9.534663309335294E-001,2.457973973823296E-002},
    {9.748463285901534E-001,1.815657770961400E-002},
    {9.897394542663850E-001,1.161344471646849E-002},
    {9.980499305356873E-001,5.002880749638552E-003},
  };
  template<>
  point_weight_pair GaussLegendre<100>::data[50]={
    {1.562898442154312E-002,3.125542345386416E-002},
    {4.687168242159179E-002,3.122488425484931E-002},
    {7.806858281343659E-002,3.116383569621074E-002},
    {1.091892035800607E-001,3.107233742756551E-002},
    {1.402031372361136E-001,3.095047885049043E-002},
    {1.710800805386032E-001,3.079837903115330E-002},
    {2.017898640957357E-001,3.061618658397945E-002},
    {2.323024818449742E-001,3.040407952645461E-002},
    {2.625881203715034E-001,3.016226510516754E-002},
    {2.926171880384720E-001,2.989097959333334E-002},
    {3.223603439005293E-001,2.959048805991241E-002},
    {3.517885263724220E-001,2.926108411063884E-002},
    {3.808729816246301E-001,2.890308960112587E-002},
    {4.095852916783023E-001,2.851685432239523E-002},
    {4.378974021720319E-001,2.810275565910164E-002},
    {4.657816497733580E-001,2.766119822079209E-002},
    {4.932107892081916E-001,2.719261344657662E-002},
    {5.201580198817631E-001,2.669745918357103E-002},
    {5.465970120650943E-001,2.617621923954466E-002},
    {5.725019326213812E-001,2.562940291020853E-002},
    {5.978474702471789E-001,2.505754448157973E-002},
    {6.226088602037082E-001,2.446120270795577E-002},
    {6.467619085141294E-001,2.384096026596899E-002},
    {6.702830156031415E-001,2.319742318525421E-002},
    {6.931491993558019E-001,2.253122025633597E-002},
    {7.153381175730568E-001,2.184300241624767E-002},
    {7.368280898020211E-001,2.113344211252747E-002},
    {7.575981185197069E-001,2.040323264620976E-002},
    {7.776279096494952E-001,1.965308749443533E-002},
    {7.968978923903148E-001,1.888373961337431E-002},
    {8.153892383391762E-001,1.809594072212793E-002},
    {8.330838798884006E-001,1.729046056832494E-002},
    {8.499645278795909E-001,1.646808617614444E-002},
    {8.660146884971646E-001,1.562962107754551E-002},
    {8.812186793850185E-001,1.477588452744186E-002},
    {8.955616449707267E-001,1.390771070371846E-002},
    {9.090295709825298E-001,1.302594789297084E-002},
    {9.216092981453337E-001,1.213145766298006E-002},
    {9.332885350430794E-001,1.122511402318562E-002},
    {9.440558701362557E-001,1.030780257486862E-002},
    {9.539007829254919E-001,9.380419653694766E-003},
    {9.628136542558153E-001,8.443871469668482E-003},
    {9.707857757637066E-001,7.499073255464270E-003},
    {9.778093584869182E-001,6.546948450845874E-003},
    {9.838775407060567E-001,5.588428003865582E-003},
    {9.889843952429919E-001,4.624450063422539E-003},
    {9.931249370374432E-001,3.655961201326225E-003},
    {9.962951347331254E-001,2.683925371553286E-003},
    {9.984919506395956E-001,1.709392653518899E-003},
    {9.997137267734412E-001,7.346344905058351E-004},
  };
  template<>
  point_weight_pair GaussLegendre<200>::data[100]={
    {7.834291142306179E-003,1.566826171583293E-002},
    {2.350095006145761E-002,1.566441506361033E-002},
    {3.916183935642233E-002,1.565672270354463E-002},
    {5.481311418495251E-002,1.564518652415393E-002},
    {7.045093206520960E-002,1.562980935763785E-002},
    {8.607145381911593E-002,1.561059497918678E-002},
    {1.016708445148974E-001,1.558754810604198E-002},
    {1.172452744085828E-001,1.556067439635545E-002},
    {1.327909198842238E-001,1.552998044777555E-002},
    {1.483039643926216E-001,1.549547379585715E-002},
    {1.637805993883118E-001,1.545716291218711E-002},
    {1.792170252645805E-001,1.541505720232031E-002},
    {1.946094522862963E-001,1.536916700345719E-002},
    {2.099541015203097E-001,1.531950358191882E-002},
    {2.252472057632053E-001,1.526607913037382E-002},
    {2.404850104661713E-001,1.520890676484301E-002},
    {2.556637746567648E-001,1.514800052148831E-002},
    {2.707797718573424E-001,1.508337535315473E-002},
    {2.858292909999357E-001,1.501504712571098E-002},
    {3.008086373373349E-001,1.494303261414732E-002},
    {3.157141333501773E-001,1.486734949845400E-002},
    {3.305421196497962E-001,1.478801635929416E-002},
    {3.452889558766267E-001,1.470505267342865E-002},
    {3.599510215939304E-001,1.461847880893864E-002},
    {3.745247171766409E-001,1.452831602022822E-002},
    {3.890064646950869E-001,1.443458644280398E-002},
    {4.033927087934002E-001,1.433731308783669E-002},
    {4.176799175623751E-001,1.423651983652022E-002},
    {4.318645834065729E-001,1.413223143420121E-002},
    {4.459432239054582E-001,1.402447348431180E-002},
    {4.599123826683536E-001,1.391327244207101E-002},
    {4.737686301830055E-001,1.379865560800386E-002},
    {4.875085646575487E-001,1.368065112122816E-002},
    {5.011288128556678E-001,1.355928795255699E-002},
    {5.146260309247439E-001,1.343459589737457E-002},
    {5.279969052167965E-001,1.330660556832780E-002},
    {5.412381531019976E-001,1.317534838781152E-002},
    {5.543465237745811E-001,1.304085658025184E-002},
    {5.673187990509334E-001,1.290316316419256E-002},
    {5.801517941596780E-001,1.276230194419404E-002},
    {5.928423585235543E-001,1.261830750253189E-002},
    {6.053873765329050E-001,1.247121519070636E-002},
    {6.177837683105758E-001,1.232106112076043E-002},
    {6.300284904680448E-001,1.216788215642190E-002},
    {6.421185368525914E-001,1.201171590404337E-002},
    {6.540509392853259E-001,1.185260070337966E-002},
    {6.658227682898954E-001,1.169057561816401E-002},
    {6.774311338116869E-001,1.152568042653329E-002},
    {6.888731859273536E-001,1.135795561124140E-002},
    {7.001461155444900E-001,1.118744234973778E-002},
    {7.112471550912781E-001,1.101418250405669E-002},
    {7.221735791959435E-001,1.083821861052681E-002},
    {7.329227053558542E-001,1.065959386934211E-002},
    {7.434918945960871E-001,1.047835213395018E-002},
    {7.538785521173128E-001,1.029453790028522E-002},
    {7.640801279328395E-001,1.010819629584683E-002},
    {7.740941174946421E-001,9.919373068619170E-003},
    {7.839180623082518E-001,9.728114575840658E-003},
    {7.935495505363276E-001,9.534467772619550E-003},
    {8.029862175907783E-001,9.338480200414964E-003},
    {8.122257467132833E-001,9.140199975351907E-003},
    {8.212658695440719E-001,8.939675776422485E-003},
    {8.301043666788187E-001,8.736956833527362E-003},
    {8.387390682135248E-001,8.532092915385720E-003},
    {8.471678542772392E-001,8.325134317331827E-003},
    {8.553886555525035E-001,8.116131848949227E-003},
    {8.633994537833832E-001,7.905136821609989E-003},
    {8.711982822709637E-001,7.692201035873163E-003},
    {8.787832263561884E-001,7.477376768768350E-003},
    {8.861524238899223E-001,7.260716760958164E-003},
    {8.933040656901267E-001,7.042274203802085E-003},
    {9.002363959860233E-001,6.822102726284827E-003},
    {9.069477128491587E-001,6.600256381861296E-003},
    {9.134363686112408E-001,6.376789635183933E-003},
    {9.197007702686635E-001,6.151757348730323E-003},
    {9.257393798736083E-001,5.925214769351110E-003},
    {9.315507149116372E-001,5.697217514688546E-003},
    {9.371333486656787E-001,5.467821559550593E-003},
    {9.424859105663231E-001,5.237083222156630E-003},
    {9.476070865283423E-001,5.005059150339589E-003},
    {9.524956192733554E-001,4.771806307640666E-003},
    {9.571503086385662E-001,4.537381959359005E-003},
    {9.615700118715106E-001,4.301843658513322E-003},
    {9.657536439107449E-001,4.065249231775649E-003},
    {9.697001776524374E-001,3.827656765344403E-003},
    {9.734086442028306E-001,3.589124590811235E-003},
    {9.768781331165662E-001,3.349711271035094E-003},
    {9.801077926209241E-001,3.109475586102519E-003},
    {9.830968298260959E-001,2.868476519469471E-003},
    {9.858445109217865E-001,2.626773244516972E-003},
    {9.883501613607615E-001,2.384425112004731E-003},
    {9.906131660306734E-001,2.141491639441064E-003},
    {9.926329694171897E-001,1.898032504938859E-003},
    {9.944090757657010E-001,1.654107552313073E-003},
    {9.959410492610110E-001,1.409776827453116E-003},
    {9.972285142833798E-001,1.165100714756534E-003},
    {9.982711559489105E-001,9.201404593416220E-004},
    {9.990687218731519E-001,6.749606344802847E-004},
    {9.996210312809368E-001,4.296466304500147E-004},
    {9.999280712850698E-001,1.845900974717014E-004},
  };

  template<>
  point_weight_pair GaussLegendre<8>::data[4]={
    {1.834346424956496E-001,3.626837833783625E-001},
    {5.255324099163290E-001,3.137066458778878E-001},
    {7.966664774136266E-001,2.223810344533740E-001},
    {9.602898564975362E-001,1.012285362903756E-001},
  };
  template<>
  point_weight_pair GaussLegendre<16>::data[8]={
    {9.501250983763747E-002,1.894506104550681E-001},
    {2.816035507792589E-001,1.826034150449247E-001},
    {4.580167776572275E-001,1.691565193950028E-001},
    {6.178762444026439E-001,1.495959888165762E-001},
    {7.554044083550030E-001,1.246289712555344E-001},
    {8.656312023878315E-001,9.515851168249333E-002},
    {9.445750230732319E-001,6.225352393864743E-002},
    {9.894009349916498E-001,2.715245941175365E-002},
  };
  template<>
  point_weight_pair GaussLegendre<32>::data[16]={
    {4.830766568773810E-002,9.654008851472772E-002},
    {1.444719615827966E-001,9.563872007927453E-002},
    {2.392873622521368E-001,9.384439908080457E-002},
    {3.318686022821275E-001,9.117387869576342E-002},
    {4.213512761306357E-001,8.765209300440363E-002},
    {5.068999089322295E-001,8.331192422694644E-002},
    {5.877157572407624E-001,7.819389578707035E-002},
    {6.630442669302150E-001,7.234579410884842E-002},
    {7.321821187402895E-001,6.582222277636177E-002},
    {7.944837959679423E-001,5.868409347853618E-002},
    {8.493676137325699E-001,5.099805926237627E-002},
    {8.963211557660517E-001,4.283589802222663E-002},
    {9.349060759377400E-001,3.427386291302144E-002},
    {9.647622555875065E-001,2.539206530926234E-002},
    {9.856115115452684E-001,1.627439473090517E-002},
    {9.972638618494815E-001,7.018610009470471E-003},
  };
  template<>
  point_weight_pair GaussLegendre<64>::data[32]={
    {2.435029266342412E-002,4.869095700913995E-002},
    {7.299312178779892E-002,4.857546744150346E-002},
    {1.214628192961206E-001,4.834476223480341E-002},
    {1.696444204239931E-001,4.799938859645788E-002},
    {2.174236437400074E-001,4.754016571483066E-002},
    {2.646871622087672E-001,4.696818281621022E-002},
    {3.113228719902110E-001,4.628479658131266E-002},
    {3.572201583376682E-001,4.549162792741684E-002},
    {4.022701579639919E-001,4.459055816375725E-002},
    {4.463660172534639E-001,4.358372452932376E-002},
    {4.894031457070531E-001,4.247351512365347E-002},
    {5.312794640198949E-001,4.126256324262433E-002},
    {5.718956462026343E-001,3.995374113272104E-002},
    {6.111553551723933E-001,3.855015317861573E-002},
    {6.489654712546571E-001,3.705512854023933E-002},
    {6.852363130542336E-001,3.547221325688256E-002},
    {7.198818501716106E-001,3.380516183714313E-002},
    {7.528199072605321E-001,3.205792835485115E-002},
    {7.839723589433414E-001,3.023465707240144E-002},
    {8.132653151227981E-001,2.833967261425893E-002},
    {8.406292962525808E-001,2.637746971505373E-002},
    {8.659993981540931E-001,2.435270256871282E-002},
    {8.893154459951143E-001,2.227017380838316E-002},
    {9.105221370785034E-001,2.013482315352975E-002},
    {9.295691721319400E-001,1.795171577569733E-002},
    {9.464113748584031E-001,1.572603047602591E-002},
    {9.610087996520538E-001,1.346304789671876E-002},
    {9.733268277899109E-001,1.116813946013003E-002},
    {9.833362538846262E-001,8.846759826363879E-003},
    {9.910133714767445E-001,6.504457968978086E-003},
    {9.963401167719554E-001,4.147033260561887E-003},
    {9.993050417357726E-001,1.783280721697018E-003},
  };
  template<>
  point_weight_pair GaussLegendre<128>::data[64]={
    {1.222369896061569E-002,2.444618019626333E-002},
    {3.666379096873352E-002,2.443156909785045E-002},
    {6.108196960413988E-002,2.440235563384914E-002},
    {8.546364050451531E-002,2.435855726468987E-002},
    {1.097942311276437E-001,2.430020016797132E-002},
    {1.340591994611876E-001,2.422731922281499E-002},
    {1.582440427142252E-001,2.413995798901864E-002},
    {1.823343059853374E-001,2.403816868102399E-002},
    {2.063155909020797E-001,2.392201213670340E-002},
    {2.301735642266601E-001,2.379155778100401E-002},
    {2.538939664226942E-001,2.364688358444940E-002},
    {2.774626201779049E-001,2.348807601653531E-002},
    {3.008654388776775E-001,2.331522999406260E-002},
    {3.240884350244135E-001,2.312844882438751E-002},
    {3.471177285976355E-001,2.292784414368663E-002},
    {3.699395553498590E-001,2.271353585023721E-002},
    {3.925402750332674E-001,2.248565203274682E-002},
    {4.149063795522746E-001,2.224432889379965E-002},
    {4.370245010371038E-001,2.198971066845983E-002},
    {4.588814198335522E-001,2.172194953805202E-002},
    {4.804640724041720E-001,2.144120553920906E-002},
    {5.017595591361440E-001,2.114764646822043E-002},
    {5.227551520511754E-001,2.084144778075006E-002},
    {5.434383024128104E-001,2.052279248696168E-002},
    {5.637966482266179E-001,2.019187104212900E-002},
    {5.838180216287635E-001,1.984888123283098E-002},
    {6.034904561585485E-001,1.949402805870801E-002},
    {6.228021939105846E-001,1.912752360995112E-002},
    {6.417416925623072E-001,1.874958694054478E-002},
    {6.602976322726462E-001,1.836044393732995E-002},
    {6.784589224477190E-001,1.796032718500812E-002},
    {6.962147083695145E-001,1.754947582711836E-002},
    {7.135543776835871E-001,1.712813542311017E-002},
    {7.304675667419087E-001,1.669655780158888E-002},
    {7.469441667970619E-001,1.625500090978559E-002},
    {7.629743300440944E-001,1.580372865939941E-002},
    {7.785484755064120E-001,1.534301076886360E-002},
    {7.936572947621934E-001,1.487312260214747E-002},
    {8.082917575079137E-001,1.439434500416662E-002},
    {8.224431169556438E-001,1.390696413295216E-002},
    {8.361029150609066E-001,1.341127128861556E-002},
    {8.492629875779687E-001,1.290756273926651E-002},
    {8.619154689395485E-001,1.239613954395104E-002},
    {8.740527969580321E-001,1.187730737273966E-002},
    {8.856677173453974E-001,1.135137632408092E-002},
    {8.967532880491582E-001,1.081866073950492E-002},
    {9.073028834017566E-001,1.027947901583308E-002},
    {9.173101980809604E-001,9.734153415007140E-003},
    {9.267692508789479E-001,9.183009871661022E-003},
    {9.356743882779164E-001,8.626377798616908E-003},
    {9.440202878302197E-001,8.064589890485766E-003},
    {9.518019613412642E-001,7.497981925634684E-003},
    {9.590147578537002E-001,6.926892566898878E-003},
    {9.656543664319653E-001,6.351663161707560E-003},
    {9.717168187471366E-001,5.772637542866070E-003},
    {9.771984914639074E-001,5.190161832675313E-003},
    {9.820961084357184E-001,4.604584256703913E-003},
    {9.864067427245860E-001,4.016254983737921E-003},
    {9.901278184917344E-001,3.425526040909731E-003},
    {9.932571129002129E-001,2.832751471458588E-003},
    {9.957927585349813E-001,2.238288430962964E-003},
    {9.977332486255137E-001,1.642503018668767E-003},
    {9.990774599773759E-001,1.045812679339914E-003},
    {9.998248879471320E-001,4.493809602915671E-004},
  };
  template<>
  point_weight_pair GaussLegendre<256>::data[128]={
    {6.123912375188530E-003,1.224767164029036E-002},
    {1.837081847881305E-002,1.224583436974840E-002},
    {3.061496877997888E-002,1.224216010427313E-002},
    {4.285452653637852E-002,1.223664939503972E-002},
    {5.508765569463373E-002,1.222930306871076E-002},
    {6.731252116571525E-002,1.222012222730384E-002},
    {7.952728910023263E-002,1.220910824803759E-002},
    {9.173012716351851E-002,1.219626278311467E-002},
    {1.039192048105096E-001,1.218158775948133E-002},
    {1.160926935603327E-001,1.216508537853558E-002},
    {1.282487672706073E-001,1.214675811579305E-002},
    {1.403856024113762E-001,1.212660872052791E-002},
    {1.525013783386553E-001,1.210464021534126E-002},
    {1.645942775675534E-001,1.208085589572390E-002},
    {1.766624860449020E-001,1.205525932956004E-002},
    {1.887041934213889E-001,1.202785435658299E-002},
    {2.007175933231267E-001,1.199864508780443E-002},
    {2.127008836226264E-001,1.196763590490547E-002},
    {2.246522667091317E-001,1.193483145956444E-002},
    {2.365699497582840E-001,1.190023667276683E-002},
    {2.484521450010568E-001,1.186385673407176E-002},
    {2.602970699919426E-001,1.182569710082494E-002},
    {2.721029478763365E-001,1.178576349734342E-002},
    {2.838680076570816E-001,1.174406191405965E-002},
    {2.955904844601360E-001,1.170059860661985E-002},
    {3.072686197993195E-001,1.165538009494373E-002},
    {3.189006618401070E-001,1.160841316225339E-002},
    {3.304848656624174E-001,1.155970485404371E-002},
    {3.420194935223719E-001,1.150926247704015E-002},
    {3.535028151129696E-001,1.145709359809150E-002},
    {3.649331078236542E-001,1.140320604303741E-002},
    {3.763086569987166E-001,1.134760789554634E-002},
    {3.876277561945155E-001,1.129030749587536E-002},
    {3.988887074354590E-001,1.123131343964957E-002},
    {4.100898214687165E-001,1.117063457655232E-002},
    {4.212294180176238E-001,1.110828000900938E-002},
    {4.323058260337412E-001,1.104425909081373E-002},
    {4.433173839475274E-001,1.097858142573018E-002},
    {4.542624399175902E-001,1.091125686604924E-002},
    {4.651393520784798E-001,1.084229551111464E-002},
    {4.759464887869835E-001,1.077170770580474E-002},
    {4.866822288668906E-001,1.069950403897982E-002},
    {4.973449618521820E-001,1.062569534189639E-002},
    {5.079330882286158E-001,1.055029268658180E-002},
    {5.184450196736749E-001,1.047330738416961E-002},
    {5.288791792948221E-001,1.039475098321283E-002},
    {5.392340018660595E-001,1.031463526793356E-002},
    {5.495079340627187E-001,1.023297225647877E-002},
    {5.596994346944810E-001,1.014977419909506E-002},
    {5.698069749365685E-001,1.006505357630777E-002},
    {5.798290385590825E-001,9.978823097035915E-003},
    {5.897641221544545E-001,9.891095696695497E-003},
    {5.996107353629682E-001,9.801884535257493E-003},
    {6.093674010963337E-001,9.711202995266292E-003},
    {6.190326557592608E-001,9.619064679840597E-003},
    {6.286050494690147E-001,9.525483410628940E-003},
    {6.380831462729114E-001,9.430473225737446E-003},
    {6.474655243637246E-001,9.334048377624339E-003},
    {6.567507762929730E-001,9.236223330955974E-003},
    {6.659375091820483E-001,9.137012760451719E-003},
    {6.750243449311628E-001,9.036431548661511E-003},
    {6.840099204260758E-001,8.934494783758384E-003},
    {6.928928877425771E-001,8.831217757249123E-003},
    {7.016719143486846E-001,8.726615961698799E-003},
    {7.103456833045436E-001,8.620705088399565E-003},
    {7.189128934599715E-001,8.513501025022425E-003},
    {7.273722596496525E-001,8.405019853220811E-003},
    {7.357225128859178E-001,8.295277846236376E-003},
    {7.439624005491116E-001,8.184291466439182E-003},
    {7.520906865754921E-001,8.072077362874030E-003},
    {7.601061516426553E-001,7.958652368754968E-003},
    {7.680075933524456E-001,7.844033498939482E-003},
    {7.757938264113259E-001,7.728237947381056E-003},
    {7.834636828081841E-001,7.611283084545282E-003},
    {7.910160119895461E-001,7.493186454806755E-003},
    {7.984496810321705E-001,7.373965773813512E-003},
    {8.057635748129983E-001,7.253638925833220E-003},
    {8.129565961764316E-001,7.132223961075179E-003},
    {8.200276660989175E-001,7.009739092970388E-003},
    {8.269757238508125E-001,6.886202695447578E-003},
    {8.337997271555051E-001,6.761633300174624E-003},
    {8.404986523457624E-001,6.636049593780484E-003},
    {8.470714945172967E-001,6.509470415052716E-003},
    {8.535172676795030E-001,6.381914752107581E-003},
    {8.598350049033763E-001,6.253401739542640E-003},
    {8.660237584665548E-001,6.123950655568467E-003},
    {8.720825999954881E-001,5.993580919114516E-003},
    {8.780106206047067E-001,5.862312086921393E-003},
    {8.838069310331580E-001,5.730163850601694E-003},
    {8.894706617776110E-001,5.597156033681580E-003},
    {8.950009632230843E-001,5.463308588644541E-003},
    {9.003970057703037E-001,5.328641593916064E-003},
    {9.056579799601445E-001,5.193175250869292E-003},
    {9.107830965950652E-001,5.056929880786590E-003},
    {9.157715868574903E-001,4.919925921814846E-003},
    {9.206227024251463E-001,4.782183925892804E-003},
    {9.253357155833165E-001,4.643724555678898E-003},
    {9.299099193340058E-001,4.504568581448246E-003},
    {9.343446275020032E-001,4.364736877968642E-003},
    {9.386391748378148E-001,4.224250421381676E-003},
    {9.427929171174623E-001,4.083130286053916E-003},
    {9.468052312391277E-001,3.941397641409080E-003},
    {9.506755153166278E-001,3.799073748766313E-003},
    {9.544031887697163E-001,3.656179958141752E-003},
    {9.579876924111783E-001,3.512737705055407E-003},
    {9.614284885307322E-001,3.368768507315766E-003},
    {9.647250609757063E-001,3.224293961794421E-003},
    {9.678769152284898E-001,3.079335741198837E-003},
    {9.708835784807429E-001,2.933915590829914E-003},
    {9.737445997043708E-001,2.788055325327707E-003},
    {9.764595497192342E-001,2.641776825427641E-003},
    {9.790280212576217E-001,2.495102034704360E-003},
    {9.814496290254646E-001,2.348052956326620E-003},
    {9.837240097603153E-001,2.200651649839813E-003},
    {9.858508222861261E-001,2.052920227965489E-003},
    {9.878297475648608E-001,1.904880853499463E-003},
    {9.896604887450652E-001,1.756555736330630E-003},
    {9.913427712075832E-001,1.607967130749461E-003},
    {9.928763426088223E-001,1.459137333310463E-003},
    {9.942609729224098E-001,1.310088681901964E-003},
    {9.954964544810965E-001,1.160843557567732E-003},
    {9.965826020233817E-001,1.011424393209112E-003},
    {9.975192527567210E-001,8.618537014203115E-004},
    {9.983062664730065E-001,7.121541634733089E-004},
    {9.989435258434087E-001,5.623489540322297E-004},
    {9.994309374662617E-001,4.124632544261141E-004},
    {9.997684374092631E-001,2.625349442963582E-004},
    {9.999560500189922E-001,1.127890178222330E-004},
  };
}

//NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
}