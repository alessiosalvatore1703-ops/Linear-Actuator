#include <Arduino.h>
int IN11 = 22;
int IN12 = 23;
int IN13 = 24;
int IN14 = 25;
int IN21 = 26;
int IN22 = 27;
int IN23 = 28;
int IN24 = 29;
int IN31 = 30;
int IN32 = 31;
int IN33 = 32;
int IN34 = 33;
int IN41 = 34;
int IN42 = 35;
int IN43 = 36;
int IN44 = 37;
int IN51 = 38;
int IN52 = 39;
int IN53 = 40;
int IN54 = 41;
int IN61 = 42;
int IN62 = 43;
int IN63 = 44;
int IN64 = 45;
int IN71 = 46;
int IN72 = 47;
int IN73 = 48;
int IN74 = 49;
int IN81 = 50;
int IN82 = 51;
int IN83 = 52;
int IN84 = 53;


struct Solenoide {
    int pin1;
    int pin2;
};

Solenoide matriceSolenoidi[2][8] = {
    // Riga 1 
    {{IN14, IN13}, {IN24, IN23}, {IN34, IN33}, {IN44, IN43}, 
     {IN54, IN53}, {IN64, IN63}, {IN74, IN73}, {IN84, IN83}},
    // Riga 2 
    {{IN11, IN12}, {IN21, IN22}, {IN31, IN32}, {IN41, IN42}, 
     {IN51, IN52}, {IN61, IN62}, {IN71, IN72}, {IN81, IN82}}
};

// Create a union to easily convert float to byte
typedef union{
  float number;
  uint8_t bytes[4];
} FLOATUNION_t;

// Create the variable you want to send
FLOATUNION_t myValue;
FLOATUNION_t xUnion;

float prec_foward = -INFINITY;
float prec_backwards = INFINITY;
float counter_foward = 0;
float counter_backwards = 0;
float dist_giri = 0;
float x = 0;
const float pitch_coils = 1;
const int num_magn = 4;
const float diam_coil = 25;
const int num_coils = 8;
const float pitch = diam_coil + pitch_coils;
unsigned long previousMillis = 0;
unsigned long loopCounter = 0;
const unsigned long LOOP_PERIOD = 10;  // Periodo del loop in ms (20Hz)

float mover[2][num_magn] = { 
    {0, 39, 72, 110},              
    {16, 39 + 16, 72 + 16, 110 + 16}            
};

const float coils[num_coils] = {
    pitch * 8, pitch * 9, pitch * 10, pitch * 11,
    pitch * 12, pitch * 13, pitch * 14, pitch * 15
};


/*const float F_ATT_ARRAY[] = {
    0.0, 0.187021081820441, 0.353159506223901, 0.500535346256711, 0.631268674965201, 0.747479565395704, 0.851288090594551, 0.944814323608073, 1.03017833748260, 1.10950020526447,
    1.18490000000000, 1.25849779473553, 1.33241366251740, 1.40876767639193, 1.48967990940545, 1.57727043460430, 1.67365932503480, 1.78096665374329, 1.90131249377610, 2.03681691817956,
    2.18960000000000, 2.36086713923742, 2.54816504370650, 2.74812574817558, 2.95738128741300, 3.17256369618711, 3.39030500926625, 3.60723726141877, 3.81999248741300, 4.02520272201730,
    4.21950000000000, 4.40039274831478, 4.56889496265660, 4.72689703090575, 4.87628934094254, 5.01896228064725, 5.15680623790019, 5.29171160058164, 5.42556875657189, 5.56026809375125,
    5.69770000000000, 5.83920076750348, 5.98389030566712, 6.13033442820141, 6.27709894881684, 6.42274968122387, 6.56585243913300, 6.70497303625469, 6.83867728629944, 6.96553100297771,
    7.08410000000000, 7.19249768167133, 7.28702781467493, 7.36354175628860, 7.41789086379011, 7.44592649445726, 7.44350000556783, 7.40646275439960, 7.33066609823036, 7.21196139433790,
    7.04620000000000, 6.83135000581121, 6.57384643563315, 6.28224104664419, 5.96508559602271, 5.63093184094709, 5.28833153859570, 4.94583644614692, 4.61199832077912, 4.29536891967069,
    4.00450000000000, 3.74649919508381, 3.52269764279245, 3.33298235713462, 3.17724035211902, 3.05535864175437, 2.96722424004937, 2.91272416101273, 2.89174541865315, 2.90417502697934,
    2.94990000000000, 3.02803861385353, 3.13463419319704, 3.26496132481733, 3.41429459550119, 3.57790859203541, 3.75107790120681, 3.92907710980217, 4.10718080460829, 4.28066357241197,
    4.44480000000000, 4.59535724950207, 4.73007278441939, 4.84717664359608, 4.94489886587624, 5.02146949010397, 5.07511855512339, 5.10407609977859, 5.10657216291369, 5.08083678337279,
    5.02510000000000, 4.93871568813820, 4.82553306912539, 4.69052520079836, 4.53866514099387, 4.37492594754871, 4.20428067829965, 4.03170239108347, 3.86216414373695, 3.70063899409687,
    3.55210000000000, 3.42043299794512, 3.30517493907903, 3.20477555321048, 3.11768457014828, 3.04235171970120, 2.97722673167802, 2.92075933588753, 2.87139926213851, 2.82759624023974,
    2.78780000000000, 2.75099312008129, 2.71828957455850, 2.69133618635970, 2.67177977841300, 2.66126717364649, 2.66144519498825, 2.67396066536640, 2.70046040770901, 2.74259124494417,
    2.80200000000000, 2.87947042172970, 2.97233396268699, 3.07705900135071, 3.19011391619972, 3.30796708571286, 3.42708688836896, 3.54394170264688, 3.65499990702547, 3.75672987998356,
    3.84560000000000, 3.91945719299992, 3.98166257469356, 4.03695580823745, 4.09007655678812, 4.14576448350209, 4.20875925153590, 4.28380052404607, 4.37562796418912, 4.48898123512159,
    4.62860000000000, 4.79766210627063, 4.99309813853878, 5.21027686569949, 5.44456705664780, 5.69133748027877, 5.94595690548743, 6.20379410116884, 6.46021783621805, 6.71059687953008,
    6.95030000000000, 7.17497398191756, 7.38137767115133, 7.56654792896461, 7.72752161662068, 7.86133559538284, 7.96502672651437, 8.03563187127856, 8.07018789093870, 8.06573164675809,
    8.01930000000000, 7.92910386605915, 7.79805037685591, 7.63022071844209, 7.42969607686949, 7.20055763818989, 6.94688658845510, 6.67276411371692, 6.38227140002714, 6.07948963343758,
    5.76850000000000, 5.45321405384585, 5.13686482142502, 4.82251569726702, 4.51323007590139, 4.21207135185761, 3.92210291966523, 3.64638817385377, 3.38799050895272, 3.14997331949163,
    2.93540000000000, 2.74637821855745, 2.58119273744402, 2.43717259248982, 2.31164681952500, 2.20194445437967, 2.10539453288396, 2.01932609086802, 1.94106816416195, 1.86794978859591,
    1.79730000000000, 1.72690777192434, 1.65640182879891, 1.58587083277369, 1.51540344599864, 1.44508833062373, 1.37501414879891, 1.30526956267417, 1.23594323439946, 1.16712382612475,
    1.09890000000000, 1.03136091374519, 0.964597707360326, 0.898702016415403, 0.833765476480435, 0.769879723125424, 0.707136391920379, 0.645627118435309, 0.585443538240217, 0.526677286905113,
    0.469420000000000, 0.413763313094887, 0.359798861759783, 0.307618281564691, 0.257313208079621, 0.208975276874576, 0.162696123519565, 0.118567383584597, 0.0766806926396740, 0.0371276862548072,
    0.0
};*/


/*const float F_REP_ARRAY[] = {
    0.0, -0.0134930950589741, -0.0169259706256407, -0.0107737844213671, 0.00448830583247912, 0.0283851424145304, 0.0604415676034193, 0.100182423677778, 0.147132552916240, 0.200816797597436,
    0.260760000000000, 0.326487002402564, 0.397522647083761, 0.473391776322222, 0.553619232396581, 0.637729857585470, 0.725248494167521, 0.815699984421367, 0.908609170625641, 1.00350089505897,
    1.09990000000000, 1.19738086544872, 1.29571602229060, 1.39472753913248, 1.49423748458120, 1.59406792724359, 1.69404093572650, 1.79397857863675, 1.89370292458120, 1.99303604216667,
    2.09180000000000, 2.18964987580256, 2.28557278375384, 2.37838884714786, 2.46691818927863, 2.54998093344017, 2.62639720292649, 2.69498712103162, 2.75457081104957, 2.80396839627436,
    2.84200000000000, 2.86803893134103, 2.88367124269403, 2.89103617227608, 2.89227295830429, 2.88952083899575, 2.88491905256754, 2.88060683723677, 2.87872343122053, 2.88140807273591,
    2.89080000000000, 2.90748939883332, 2.92587024547005, 2.93878746374781, 2.93908597750421, 2.91961071057684, 2.87320658680333, 2.79271853002129, 2.67099146406831, 2.50087031278201,
    2.27520000000000, 1.99022907332570, 1.65582057542576, 1.28524117273267, 0.891757531678884, 0.488636318696879, 0.0891442002191205, -0.293452157321915, -0.645886087493763, -0.954890923863947,
    -1.20720000000000, -1.39277929213611, -1.51452534717312, -1.57856735467850, -1.59103450421975, -1.55805598536436, -1.48576098767982, -1.38027870073363, -1.24773831409327, -1.09426901732623,
    -0.926000000000000, -0.748495004781248, -0.565055986733305, -0.378419454018685, -0.191321914799894, -0.00649987723943968, 0.173310150500170, 0.345371660256430, 0.506948143866825, 0.655303093168851,
    0.787700000000000, 0.901899111261102, 0.997647694106343, 1.07518977075324, 1.13476936341933, 1.17663049432212, 1.20101718567914, 1.20817345970792, 1.19834333862597, 1.17177084465082,
    1.12870000000000, 1.06985573973684, 0.997886650307939, 0.915922231005718, 0.827091981122589, 0.734525399950967, 0.641351986783266, 0.550701240911905, 0.465702661629297, 0.389485748227857,
    0.325180000000000, 0.275177409791541, 0.238919944661898, 0.215112065223885, 0.202458232090318, 0.199662905874015, 0.205430547187790, 0.218465616644460, 0.237472574856840, 0.261155882437748,
    0.288220000000000, 0.317908691096994, 0.351622931044468, 0.391302998098743, 0.438889170516138, 0.496321726552975, 0.565540944465574, 0.648487102510255, 0.747100478943338, 0.863321352021146,
    0.999090000000000, 1.15531533582048, 1.32878081116023, 1.51523851238114, 1.71044052584513, 1.91013893791409, 2.11008583494991, 2.30603330331452, 2.49373342936980, 2.66893829947766,
    2.82740000000000, 2.96643211562108, 3.08959422431462, 3.20200740237668, 3.30879272610335, 3.41507127179068, 3.52596411573476, 3.64659233423166, 3.78207700357745, 3.93753920006821,
    4.11810000000000, 4.32726701169519, 4.56209397158129, 4.81802114811213, 5.09048880974148, 5.37493722492318, 5.66680666211102, 5.96153738975882, 6.25456967632038, 6.54134379024950,
    6.81730000000000, 7.07808923759817, 7.32020508936020, 7.54035180517483, 7.73523363493072, 7.90155482851660, 8.03601963582114, 8.13533230673305, 8.19619709114104, 8.21531823893379,
    8.18940000000000, 8.11631783791216, 7.99863207097789, 7.84007423118857, 7.64437585053563, 7.41526846101044, 7.15648359460442, 6.87175278330897, 6.56480755911548, 6.23937945401536,
    5.89920000000000, 5.54819411075319, 5.19106022672826, 4.83269017007088, 4.47797576292679, 4.13180882744165, 3.79908118576118, 3.48468466003108, 3.19351107239705, 2.93045224500480,
    2.70040000000000, 2.50688461907506, 2.34799022210909, 2.22043938852789, 2.12095469775725, 2.04625872922297, 1.99307406235086, 1.95812327656671, 1.93812895129631, 1.92981366596548,
    1.92990000000000, 1.93523601294656, 1.94317168483538, 1.95118247581757, 1.95674384604423, 1.95733125566646, 1.95042016483538, 1.93348603370209, 1.90400432241769, 1.85945049113329,
    1.79730000000000, 1.71590690913870, 1.61713967854939, 1.50374536820184, 1.37847103806585, 1.24406374811118, 1.10327055830761, 0.958838528624940, 0.813514719032923, 0.670046189501353,
    0.531180000000000, 0.399663210498647, 0.278242880967077, 0.169666071375061, 0.0766798416923846, 0.00203125188882070, 0.0515326380658482, 0.0812647682018426, 0.0844180785493857, 0.0582455091386984,
    0.0
};

// Funzione di interpolazione lineare per la repulsione
float F_att(float x_mm) {
    // Convert mm to array index (0.1mm steps)
    float x = x_mm * 10.0; // convert mm to 0.1mm units
    if(x <= 0.0 || x >= 229.0) return 0.0; // out of range
    
    int index = (int)x;
    float frac = x - index;
    
    // Linear interpolation
    return F_ATT_ARRAY[index] + (F_ATT_ARRAY[index+1] - F_ATT_ARRAY[index]) * frac;
}

float F_rep(float x_mm) {
    // Convert mm to array index (0.1mm steps)
    float x = x_mm * 10.0; // convert mm to 0.1mm units
    if(x <= 0.0 || x >= 229.0) return 0.0; // out of range
    
    int index = (int)x;
    float frac = x - index;
    
    // Linear interpolation
    return F_REP_ARRAY[index] + (F_REP_ARRAY[index+1] - F_REP_ARRAY[index]) * frac;
}*/



void setup() {
  // initialize serial, use the same boudrate in the Simulink Config block
  Serial.begin(115200);
  for (int i = 0; i < 8; i++) {
        pinMode(matriceSolenoidi[0][i].pin1, OUTPUT);
        pinMode(matriceSolenoidi[0][i].pin2, OUTPUT);
        pinMode(matriceSolenoidi[1][i].pin1, OUTPUT);
        pinMode(matriceSolenoidi[1][i].pin2, OUTPUT);
    }
}
float getFloat(){
    int cont = 0;
    FLOATUNION_t f;
    while (cont < 4 ){
        f.bytes[cont] = Serial.read() ;
        cont = cont + 1;
    }
    return f.number;
  }

// IN TEORIA ELIMINABILE 
/*int calculateCoilMode(float coil_pos, float mover_row[], float x) {
    for (int j = 0; j < num_magn; j++) {
        float dist = coil_pos - (mover_row[j] + x);
        float abs_dist = fabs(dist);
        if (abs_dist < 23) {
            int mode = (dist > 0) ? 1 : ((dist < 0) ? -1 : 0);
            if (mode != 0 && j < num_magn - 1) {
                float dist_next = coil_pos - (mover_row[j + 1] + x);
                float abs_dist_next = fabs(dist_next);
                
                if (abs_dist_next < 23) {
                    // Scegli il magnete più vicino
                    bool next_is_closer = (abs_dist_next < abs_dist);
                    float closest_pos = mover_row[next_is_closer ? j + 1 : j];
                    mode = (coil_pos - closest_pos > 0) ? 1 : -1;
                    return mode;
                }
            }
            return mode;
        }
    }
    return 0;
  }*/


/*ELIMINABILE, SFRUTTA I DATI DELLE CFORZE CHE SONO IMPRECISI
int calculateCoilMode(float coil_pos, float mover_row[], float x) {
  float F_total_rep = 0;
  float F_total_att = 0;
  
  for (int j = 0; j < num_magn; j++) {
    float dist = coil_pos - (mover_row[j] + x);
    float abs_dist = fabs(dist);
    
    if (abs_dist < 23) {  // Soglia di 23 unità
      if (dist < 0) {      // Controllo direzione
        F_total_rep -= F_rep(abs_dist);
        F_total_att -= F_att(abs_dist);
      } else {
        F_total_rep += F_rep(abs_dist);
        F_total_att += F_att(abs_dist);
      }
    }
  }
  
  // Condizione aggiunta per forze nulle
  if (F_total_att != 0 || F_total_rep != 0) {
    if (F_total_att > F_total_rep) {
      return 1;   // Attrazione
    } else {
      return -1;  // Repulsione
    }
  } else {
    return 0;     // Nessuna forza significativa
  }
}*/
//accendo e spengo solenoide in base a dove sta il magnete piu vicino
int calculateCoilMode(float coil_pos, float mover_row[], float x) {
    float min_abs_dist = 100000.0;  // grande valore iniziale
    int closest_index = -1;

    for (int j = 0; j < num_magn; j++) {
        float dist = coil_pos - (mover_row[j] + x);
        float abs_dist = fabs(dist);

        if (abs_dist < 23 && abs_dist < min_abs_dist) { //nota si puo giocare sul valore in zero, che comunque influisce per i prossimi 5 milli
            min_abs_dist = abs_dist;
            closest_index = j;
        }
    }

    if (closest_index != -1) {
        float closest_pos = mover_row[closest_index] + x;
        float dist = coil_pos - closest_pos;
        if (dist > 0) return 1;
        else if (dist < 0) return -1;
        else if (dist == 0) return 0;
    }

    return 0;
} 

void controlSolenoids() {
    for (int m = 0; m < 2; m++) {  // Per ogni riga di solenoidi
        for (int s = 0; s < num_coils; s++) {  // Per ogni solenoide
            int mode = calculateCoilMode(coils[s], mover[m], x);
            // Controllo del solenoide
            digitalWrite(matriceSolenoidi[m][s].pin1, (mode == -1) ? HIGH : LOW); //setta in maniera corretta pin 1 e 2 in base a mode
            digitalWrite(matriceSolenoidi[m][s].pin2, (mode == 1) ? HIGH : LOW);
        }
    }
  }
void loop(){
  unsigned long currentMillis = millis();
  unsigned long elapsed = currentMillis - previousMillis;
// Esegui il codice solo se è passato il periodo desiderato
    if (elapsed >= LOOP_PERIOD) {
        previousMillis = currentMillis;
        loopCounter++;
  myValue.number = getFloat(); // Acquisisci il valore dal sensore
  int alpha = 20919; // per settare in maniera corretta alpha --> quando si accende per la prima volta il sensore leggere che valore returna se posizionato in zero e inserirlo qua
  // === Calcolo posizione continua ===
  float x1 = (myValue.number  - alpha) / 8192.0; 
   /* if (x1 - prec_foward < -7.5) counter_foward++; //logica giri momentaneamente eliminata
    if (prec_backwards - x1 < -7.5) counter_backwards++; //7.5 perche ho un salto di circa 8 giri quando fa overflow
  prec_foward = x1;
  prec_backwards = x1;
  float dist_giri = (counter_foward - counter_backwards) * 48 * PI * 8;*/
  x = x1 * 48 * PI;

for (int m = 0; m < 2; m++) {  // Per ogni riga di solenoidi
  for (int s = 0; s < num_coils; s++) {  // Per ogni solenoide
            int mode = 1;
            // Controllo del solenoide
            digitalWrite(matriceSolenoidi[m][s].pin1, (mode == -1) ? HIGH : LOW); //setta in maniera corretta pin 1 e 2 in base a mode
            digitalWrite(matriceSolenoidi[m][s].pin2, (mode == 1) ? HIGH : LOW);
        }
    }




// Invio del valore x sulla seriale
    xUnion.number = x; // Carichiamo il valore x nell'union
    controlSolenoids();
    //parte eliminabile che rimanda il valore a simulink
    Serial.write('A');                 // Header
    for (int i = 0; i < 4; i++) {      // Payload
      Serial.write(xUnion.bytes[i]);    // Dato o zero
    }
    Serial.print('\n');                // Terminatore
    }
    else {
        delayMicroseconds((LOOP_PERIOD - elapsed) * 1000);
   }
//

}
