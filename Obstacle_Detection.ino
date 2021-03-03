#include <Servo.h>

// constants won't change
const int TRIG_PIN  = 10;  // Arduino pin connected to Ultrasonic Sensor's TRIG pin
const int ECHO_PIN  = 11;  // Arduino pin connected to Ultrasonic Sensor's ECHO pin
const int SERVO_PIN = 12;  // Arduino pin connected to Servo Motor's pin
const int BUZZER    = 9;   // Arduino pin connected to piezoelectric buzzer
const int DISTANCE_THRESHOLD = 50; // centimeters

Servo servo; // create servo object to control a servo

// variables will change:
float duration_us, distance_cm;

void setup() {
  Serial.begin (9600);       // initialize serial port
  pinMode(TRIG_PIN, OUTPUT); // set arduino pin to output mode
  pinMode(ECHO_PIN, INPUT);  // set arduino pin to input mode
  pinMode(BUZZER, OUTPUT);   // set buzzer pin to output mode
  servo.attach(SERVO_PIN);   // attaches the servo on pin 9 to the servo object
}

void loop() {
  Serial.println("Starting obstacle detection...");
  int direction = 0;
  // Yaw angle to rotate from accelerometer data
  int angle = 44;
  // The specs of this servo state (from datasheet):
  //          0.12s/60 degrees (4.8V)
  //          0.10s/60 degrees (6V)
  // Linear approximation to find the speed of our specific servo
  // in s/60 degrees.
  float T = (0.12 - 0.10) * (5.0 - 4.8) / (4.8 - 6.0) + 0.12;
  int turn_time = T * (float(angle) / 60.0) * 1000;
  if (angle > 0){
    servo.write(direction); // start rotating clockwise since angle is positive 
    delay(turn_time); // rotate until angle is reached
    servo.write(90); // stop rotating
  }
  else if (angle < 0){
    direction = 180;
    servo.write(direction); // start rotating counterclockwise since angle is negative
    delay(turn_time); // rotate until angle is reached
    servo.write(90); // stop rotating
  }
  else{
    servo.write(90); // don't rotate since individual is moving straight
  }

  delay(100);

  float distance = calculateDistance();
    
  if(distance_cm <= DISTANCE_THRESHOLD){
    Serial.println("Object in path!");
    digitalWrite(BUZZER, HIGH);
    delay(5000);
  }
  else {
    digitalWrite(BUZZER, LOW);
  }

  delay(100);
  
  // Reset servo position
  if (direction == 0){
    servo.write(180);
    delay(turn_time);
    servo.write(90);
  }
  else if (direction == 180){
    servo.write(0);
    delay(turn_time);
    servo.write(90);
  }

  Serial.println("Waiting before we start the loop again");
  delay(5000);
}

// Function for calculating the distance measured by the Ultrasonic sensor
int calculateDistance(){ 
  // generate 10-microsecond pulse to TRIG pin
  digitalWrite(TRIG_PIN, HIGH);
  delayMicroseconds(10);
  digitalWrite(TRIG_PIN, LOW);
  // measure duration of pulse from ECHO pin
  duration_us = pulseIn(ECHO_PIN, HIGH);
  // calculate the distance
  distance_cm = 0.017 * duration_us;
  // print the value to Serial Monitor
  Serial.print("distance: ");
  Serial.print(distance_cm);
  Serial.println(" cm");
  return distance_cm;
}
