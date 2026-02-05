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
    //config 1
    //riga1
    digitalWrite(riga1[0].pin2, LOW); 
    digitalWrite(riga1[0].pin1, HIGH);
    digitalWrite(riga1[1].pin2, LOW); 
    digitalWrite(riga1[1].pin1, HIGH);
    digitalWrite(riga1[2].pin2, HIGH); 
    digitalWrite(riga1[2].pin1, LOW);
    digitalWrite(riga1[3].pin2, HIGH); 
    digitalWrite(riga1[3].pin1, LOW);
    digitalWrite(riga1[5].pin2, HIGH); 
    digitalWrite(riga1[5].pin1, LOW);
    
    //riga2
    digitalWrite(riga2[0].pin2, LOW); 
    digitalWrite(riga2[0].pin1, HIGH);
    digitalWrite(riga2[1].pin2, LOW); 
    digitalWrite(riga2[1].pin1, HIGH);
    digitalWrite(riga2[2].pin2, HIGH); 
    digitalWrite(riga2[2].pin1, LOW); 
    digitalWrite(riga2[3].pin2, HIGH); 
    digitalWrite(riga2[3].pin1, LOW);
    digitalWrite(riga2[5].pin2, HIGH); 
    digitalWrite(riga2[5].pin1, LOW); 

    delay(38);

    //config 2
    //riga1
    digitalWrite(riga1[0].pin2, LOW); 
    digitalWrite(riga1[0].pin1, HIGH);
    digitalWrite(riga1[5].pin2, LOW); 
    digitalWrite(riga1[5].pin1, HIGH);
    digitalWrite(riga1[3].pin2, LOW); 
    digitalWrite(riga1[3].pin1, HIGH);
    digitalWrite(riga1[1].pin2, HIGH); 
    digitalWrite(riga1[1].pin1, LOW);
    digitalWrite(riga1[2].pin2, HIGH); 
    digitalWrite(riga1[2].pin1, LOW);
    digitalWrite(riga1[4].pin2, HIGH); 
    digitalWrite(riga1[4].pin1, LOW);
    
    //riga2
    digitalWrite(riga2[0].pin2, LOW); 
    digitalWrite(riga2[0].pin1, HIGH);
    digitalWrite(riga2[3].pin2, LOW); 
    digitalWrite(riga2[3].pin1, HIGH);
    digitalWrite(riga2[5].pin2, LOW); 
    digitalWrite(riga2[5].pin1, HIGH);

    digitalWrite(riga2[2].pin2, HIGH); 
    digitalWrite(riga2[2].pin1, LOW); 
    digitalWrite(riga2[1].pin2, HIGH); 
    digitalWrite(riga2[1].pin1, LOW);
    digitalWrite(riga2[4].pin2, HIGH); 
    digitalWrite(riga2[4].pin1, LOW); 

    delay(30);

    //config 3
    //riga1
    digitalWrite(riga1[2].pin2, LOW); 
    digitalWrite(riga1[2].pin1, HIGH);
    digitalWrite(riga1[5].pin2, LOW); 
    digitalWrite(riga1[5].pin1, HIGH);
    digitalWrite(riga1[1].pin2, HIGH); 
    digitalWrite(riga1[1].pin1, LOW);
    digitalWrite(riga1[4].pin2, HIGH); 
    digitalWrite(riga1[4].pin1, LOW);
    
    
    //riga2
    digitalWrite(riga2[2].pin2, LOW); 
    digitalWrite(riga2[2].pin1, HIGH);
    digitalWrite(riga2[3].pin2, LOW); 
    digitalWrite(riga2[3].pin1, HIGH);
    digitalWrite(riga2[5].pin2, LOW); 
    digitalWrite(riga2[5].pin1, HIGH);
    //digitalWrite(riga2[6].pin2, HIGH); 
    //digitalWrite(riga2[6].pin1, LOW); 
    digitalWrite(riga2[1].pin2, HIGH); 
    digitalWrite(riga2[1].pin1, LOW);
    digitalWrite(riga2[4].pin2, HIGH); 
    digitalWrite(riga2[4].pin1, LOW); 
}

