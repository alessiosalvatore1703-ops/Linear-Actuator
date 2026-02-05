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
const unsigned long LOOP_PERIOD = 50;  // Periodo del loop in ms (20Hz)

float mover[2][num_magn] = { 
    {0, 39, 72, 110},              
    {16, 39 + 16, 72 + 16, 110 + 16}            
};

const float coils[num_coils] = {
    pitch * 0, pitch * 1, pitch * 2, pitch * 3,
    pitch * 4, pitch * 5, pitch * 6, pitch * 7
};



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
   // if (elapsed >= LOOP_PERIOD) {
        previousMillis = currentMillis;
        loopCounter++;
  myValue.number = getFloat(); // Acquisisci il valore dal sensore


   // per settare in maniera corretta alpha --> quando si accende per la prima volta il sensore leggere che valore returna se posizionato in zero e inserirlo qua
  // === Calcolo posizione continua ===
  float x1 = (myValue.number) / 8192.0; 
   /* if (x1 - prec_foward < -7.5) counter_foward++; //logica giri momentaneamente eliminata
    if (prec_backwards - x1 < -7.5) counter_backwards++; //7.5 perche ho un salto di circa 8 giri quando fa overflow
  prec_foward = x1;
  prec_backwards = x1;
  float dist_giri = (counter_foward - counter_backwards) * 48 * PI * 8;*/
  x = x1 * 48 * PI;

// Invio del valore x sulla seriale
    xUnion.number = x; // Carichiamo il valore x nell'union
    controlSolenoids();
    //parte eliminabile che rimanda il valore a simulink
    Serial.write('A');                 // Header
    for (int i = 0; i < 4; i++) {      // Payload
      Serial.write(xUnion.bytes[i]);    // Dato o zero
    }
    Serial.print('\n');                // Terminatore
   // }
    //else {
       // delayMicroseconds((LOOP_PERIOD - elapsed) * 1000);
   //}
//

}
