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
int IN51 = 46;
int IN52 = 47;
int IN53 = 48;
int IN54 = 49;
int IN61 = 50;
int IN62 = 51;
int IN63 = 52;
int IN64 = 53;
struct Solenoide {
    int pin1;
    int pin2;
};

Solenoide riga1[6] = {
    {IN14, IN13},
    {IN24, IN23},
    {IN34, IN33},
    {IN44, IN43},
    {IN54, IN53},
    {IN64, IN63}
};

Solenoide riga2[6] = {
    {IN11, IN12},
    {IN21, IN22},
    {IN31, IN32},
    {IN41, IN42},
    {IN51, IN52},
    {IN61, IN62}
};

void setup() {
    for (int i = 0; i < 6; i++) {
        pinMode(riga1[i].pin1, OUTPUT);
        pinMode(riga1[i].pin2, OUTPUT);
        pinMode(riga2[i].pin1, OUTPUT);
        pinMode(riga2[i].pin2, OUTPUT);
    }
}

void loop() {
    
    // Riga 1
      unsigned long t = micros();

    if (t >= 0 && t < 38854) { digitalWrite(riga1[0].pin1, HIGH); digitalWrite(riga1[0].pin2, LOW); }
    if (t >= 0 && t < 28033) { digitalWrite(riga1[1].pin1, HIGH); digitalWrite(riga1[1].pin2, LOW); }
    if (t >= 29183 && t < 55036) { digitalWrite(riga1[1].pin1, LOW); digitalWrite(riga1[1].pin2, HIGH); }
    if (t >= 55626 && t < 67797) { digitalWrite(riga1[1].pin1, HIGH); digitalWrite(riga1[1].pin2, LOW); }
    if (t >= 0 && t < 45885) { digitalWrite(riga1[2].pin1, LOW); digitalWrite(riga1[2].pin2, HIGH); }
    if (t >= 46595 && t < 64846) { digitalWrite(riga1[2].pin1, HIGH); digitalWrite(riga1[2].pin2, LOW); }
    if (t >= 65347 && t < 67797) { digitalWrite(riga1[2].pin1, LOW); digitalWrite(riga1[2].pin2, HIGH); }
    if (t >= 0 && t < 19842) { digitalWrite(riga1[3].pin1, LOW); digitalWrite(riga1[3].pin2, HIGH); }
    if (t >= 21432 && t < 47995) { digitalWrite(riga1[3].pin1, HIGH); digitalWrite(riga1[3].pin2, LOW); }
    if (t >= 54435 && t < 67797) { digitalWrite(riga1[3].pin1, LOW); digitalWrite(riga1[3].pin2, HIGH); }
    if (t >= 35304 && t < 61756) { digitalWrite(riga1[4].pin1, LOW); digitalWrite(riga1[4].pin2, HIGH); }
    if (t >= 62276 && t < 67797) { digitalWrite(riga1[4].pin1, HIGH); digitalWrite(riga1[4].pin2, LOW); }
    if (t >= 0 && t < 36224) { digitalWrite(riga1[5].pin1, LOW); digitalWrite(riga1[5].pin2, HIGH); }
    if (t >= 37114 && t < 56796) { digitalWrite(riga1[5].pin1, HIGH); digitalWrite(riga1[5].pin2, LOW); }
    
    // Riga 2
    if (t >= 0 && t < 37924) { digitalWrite(riga2[0].pin1, HIGH); digitalWrite(riga2[0].pin2, LOW); }
    if (t >= 0 && t < 27143) { digitalWrite(riga2[1].pin1, HIGH); digitalWrite(riga2[1].pin2, LOW); }
    if (t >= 28283 && t < 53995) { digitalWrite(riga2[1].pin1, LOW); digitalWrite(riga2[1].pin2, HIGH); }
    if (t >= 54585 && t < 66937) { digitalWrite(riga2[1].pin1, HIGH); digitalWrite(riga2[1].pin2, LOW); }
    if (t >= 0 && t < 44914) { digitalWrite(riga2[2].pin1, LOW); digitalWrite(riga2[2].pin2, HIGH); }
    if (t >= 45625 && t < 63926) { digitalWrite(riga2[2].pin1, HIGH); digitalWrite(riga2[2].pin2, LOW); }
    if (t >= 64436 && t < 66937) { digitalWrite(riga2[2].pin1, LOW); digitalWrite(riga2[2].pin2, HIGH); }
    if (t >= 0 && t < 18982) { digitalWrite(riga2[3].pin1, LOW); digitalWrite(riga2[3].pin2, HIGH); }
    if (t >= 20552 && t < 47005) { digitalWrite(riga2[3].pin1, HIGH); digitalWrite(riga2[3].pin2, LOW); }
    if (t >= 53395 && t < 66937) { digitalWrite(riga2[3].pin1, LOW); digitalWrite(riga2[3].pin2, HIGH); }
    if (t >= 34383 && t < 60786) { digitalWrite(riga2[4].pin1, LOW); digitalWrite(riga2[4].pin2, HIGH); }
    if (t >= 61316 && t < 66937) { digitalWrite(riga2[4].pin1, HIGH); digitalWrite(riga2[4].pin2, LOW); }
    if (t >= 0 && t < 35294) { digitalWrite(riga2[5].pin1, LOW); digitalWrite(riga2[5].pin2, HIGH); }
    if (t >= 36194 && t < 55756) { digitalWrite(riga2[5].pin1, HIGH); digitalWrite(riga2[5].pin2, LOW); }
    
}
