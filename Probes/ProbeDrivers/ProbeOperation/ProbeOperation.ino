/*
  R – start recording.
  S – stop recording.
  P – print the recorded file. (S – stop printing the file).
  D – delete the file. (Y – to confirm deletion, N – to cancel deletion).
  F – show the files in SD/DATA/ folder.
*/

#include <SPI.h> // For communication with the SPI devices, such as SD card.
#include <SD.h> // For writing and reading files on SD card.

// Define general variables
char ch; // Character for initial Serial Monitor communication.
int numberOfFiles = 0; // Variable which keeps the total number of files in DATA folder, so that we know which file number should be next.
boolean recording = false; // A switch which tells us whether we are recoding data at the moment.
boolean showingData = false; // A switch which tells us whether we are displaying data at the moment.
boolean deletedCurrentFile = false; // A switch which detects whether the current file has been deleted.
char ansCh; //  Character for further Serial Monitor communication
int newFileNumber; // File number for recording data.
File root; // Address which tells where to write the data and where to read it from.
const int numBits = 12;

// Define variables for On/Off button operation
const int buttonPin = 2; // number for the pushbutton input pin
int buttonState = 0; // variable for reading the current pushbutton status
int lastButtonState = 0; // variable for previous pushbutton status

void setup() { // Setup code to run once only
  pinMode(buttonPin, INPUT); // initialise the pushbutton's pin to receive input

  Serial.begin(9600); // establishing connection with the serial monitor
  while (!Serial) {
    ; // wait for serial port to connect. Needed for native USB port only
  }

  Serial.print("Initializing SD card...");

  // In the initial setup pin 53 was chosen to be the communication pin.

  if (!SD.begin(53)) {
    Serial.println("initialization failed!"); //Establishing connection with SD card.
    while (1);
  }
  Serial.println("initialization done.");

  root = SD.open("/DATA"); // Declaring the path to the files.
  numberOfFiles = countFiles(root); // Counting the number of files recorded during previous sessions.
  newFileNumber = countFiles(root); // Initialising the number of the newFile. It will be increased by 1 every time user presses 'r'.
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void loop() {
  buttonState = digitalRead(buttonPin);
  // if button pressed down, and recording false
  if ((buttonState != lastButtonState) && (buttonState == HIGH) && (recording == false)) {
    lastButtonState = buttonState;
    
    Serial.println("Starting to record...");
    newFileNumber++; // Increase the number of files by one.
    deletedCurrentFile = false; // Changing the deletedCurrentFile switch, so that there is an opportunity to delete the recorded file later.
    File newFile = SD.open("/DATA/D" + String(newFileNumber) + ".txt", FILE_WRITE); // Opening file for recording data.
    analogReadResolution(numBits); // Allowing analog inputs to use 12 bits to store analog input, instead of 10 bits default.
    newFile.println("---"); // Read off the potentiometer values and save to SD file.

    Serial.println("Recording the data.");
    recording = true; // Changing the recording switch into recording mode.
    buttonState = digitalRead(buttonPin);

    while (recording) { // While we are recording…
      if (newFile) { // If the new file was opened successfully…
        newFile.print(analogRead(A0) / pow(2, numBits) * 3.3, 6); // Read off the potentiometer values and save to SD file.

        // Repeat for probes (copy and paste below 3 lines as required with A# changed on each probe added)
        delay(2); // Put a 2ms delay to increase the accuracy of readings.
        newFile.print(","); // Print the tab separation for the file to me more readable.
        newFile.print(analogRead(A1) / pow(2, numBits) * 3.3, 6); // Read off first probe value and save to SD file (6 decimal places)

        delay(2); // Put a 2ms delay to increase the accuracy of readings.
        newFile.print(","); // Print the tab separation for the file to me more readable.
        newFile.print(analogRead(A2) / pow(2, numBits) * 3.3, 6); // Read off first probe value and save to SD file (6 decimal places)
        delay(2); // Put a 2ms delay to increase the accuracy of readings.

        // Finish off each line ready for next iteration
        newFile.println(","); // Print tab and new line seperator for SD file


        // Check Button Status again
        buttonState = digitalRead(buttonPin);
        // Check for two scenarios where button status changed (other scenarios represent no change, and continue the loop)
        if ((buttonState != lastButtonState) && (buttonState == LOW)) { // When button is un-pressed, switch lastButtonState and continue loop
          lastButtonState = buttonState;
        }
        else if ((buttonState != lastButtonState) && (buttonState == HIGH)) { // when button is pressed after being un-pressed, stop recording
          recording = false; // change recording switch
          newFile.close(); // close file with recorded data
          lastButtonState = buttonState; // reset button state
          Serial.println("Done Recording"); // display end message
          File newFile = SD.open("/DATA/D" + String(newFileNumber) + ".txt"); // Open the file to read data from.
          if (newFile) { // If file opened successfully…
            while (newFile.available()) { // Run a loop until all strings are displayed or the showingData switch is changed.
              Serial.write(newFile.read()); // Display the data on the serial monitor.
              ansCh = Serial.read(); // Reading Serial Monitor input.
            }
            Serial.print("Done.");
            newFile.close(); // Close the file.
          }

        }
      }
      else { // if new file couldn't open, display error
        Serial.println("Error opening file.");
      }
    }
    delay(500); // Delay .5 seconds for stability
    // Reset variables
    buttonState = 0;
    lastButtonState = buttonState;
    recording = false;
    delay(500);
  }

  //
  // Look at other functions requiring Serial Monitor Input
  //

  if (Serial.available()) {
    ch = Serial.read(); // Read the serial monitor input
    //Case P - printing the file that has been just recorded.
    if (ch == 'P' || ch == 'p') {
      if (newFileNumber > numberOfFiles) { // If there are any files that were recorded during this session…
        if (!deletedCurrentFile) { // If the file has not been deleted…
          showingData = true; // Change the showingData switch.
          Serial.println("Showing the data that has been just recorded. ");
          File newFile = SD.open("/DATA/D" + String(newFileNumber) + ".txt"); // Open the file to read data from.
          if (newFile) { // If file opened successfully…
            while (newFile.available() && showingData) { // Run a loop until all strings are displayed or the showingData switch is changed.
              Serial.write(newFile.read()); // Display the data on the serial monitor.
              ansCh = Serial.read(); // Reading Serial Monitor input.
              if (ansCh == 's' || ansCh == 'S') {// If user uses S to stop the displaying process…
                showingData = false; // Set the displaying switch to false.
                Serial.print("\nStopped showing the data.\n");
              }
            }
            Serial.print("Done.");
            newFile.close(); // Close the file.
          }
          else Serial.println("Error opening file.");
        }
        else Serial.println("Cannot print the file. You have deleted D" + String(newFileNumber + 1) + ".txt."); // If deletedCurrentFile is true, print an error message.
      }
      else Serial.println("Cannot print the file. No files have been recorded during this session."); // If no files were recorded in this session, print an error message.
    }

    //Case D - deleting the file that has been just recorded.
    if (ch == 'd' || ch == 'D') {
      if (newFileNumber > numberOfFiles) { // If there are any files that were recorded during this session…
        if (!deletedCurrentFile) { // If the file has not been deleted…
          Serial.println("Are you sure you want to delete D" + String(newFileNumber) + ".txt? Type Y – “yes”, N – “no”."); // Display a confirmation message.
          while (true) { // Wait to get a respond.
            ansCh = Serial.read(); // Read the serial monitor input.
            delay(100); // Put delay to avoid unnecessary load of the processor.
            if (ansCh == 'y' || ansCh == 'Y') { // If Y – yes.
              SD.remove("/DATA/D" + String(newFileNumber) + ".txt"); // Delete the file.
              deletedCurrentFile = true; // Change deletedCurrentFile switch.
              Serial.println("D" + String(newFileNumber) + ".txt has been deleted successfully.");
              newFileNumber--; // Decrease the number of existing files.
              break;
            }
            if (ansCh == 'n' || ansCh == 'N') { // If N – no.
              Serial.println("Deletion of D" + String(newFileNumber) + ".txt has been cancelled."); // Print the message.
              break;
            }
          }
        }
        else Serial.println("You have already deleted D" + String(newFileNumber + 1) + ".txt."); // If deletedCurrentFile is true, print an error message.
      }
      else Serial.println("Cannot delete the file. No files have been recorded during this session."); // If no files were recorded in this session, print an error message.
    }

    //Case F - files on the sd card/data/
    if (ch == 'f' || ch == 'F') {
      Serial.println("Listing all of the files on SD/DATA/.");
      root = SD.open("/DATA/"); // Assign the address of the folder.
      printDirectory(root, 0); // Print all of the files.
      Serial.println("Done.");
    }
  }
}

// Function for printing files. Adopted from SD examples listfiles.
void printDirectory(File dir, int numTabs) {
  while (true) {
    File entry =  dir.openNextFile(); // Open a new file.
    if (! entry) {
      // No more files.
      break;
    }
    for (uint8_t i = 0; i < numTabs; i++) {
      Serial.print('\t'); // Print the tabs for better visual representation.
    }
    Serial.print(entry.name()); // Print the name of the file.
    if (entry.isDirectory()) {  // If it is a directory…
      Serial.println("/");
      printDirectory(entry, numTabs + 1); // Call the function again.
    } else {
      // Files have sizes, directories do not.
      Serial.print("\t\t");
      Serial.println(entry.size(), DEC); // Print the size of the file.
    }
    entry.close(); // Close the entry.
  }
}

//Function for counting the number of files in sd/data/.
int countFiles(File dir) {
  while (true) {
    File entry =  dir.openNextFile();
    if (! entry) {
      break;
    }
    String fileName = entry.name(); // Store the name of the file.
    if (fileName.indexOf('~') == -1) { // Check whether this is a hidden file (they normally have ~).
      numberOfFiles++; // If it is not a hidden file, add 1 to a number of files.
    }
    entry.close(); // Close the file.
  }
  return numberOfFiles; // Return the total number of files.
}
