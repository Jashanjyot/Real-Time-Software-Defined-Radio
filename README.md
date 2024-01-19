# Real-Time Software Defined Radio

This project involves the development of a real-time Software Defined Radio (SDR) using a combination of C++ and Python. The system utilizes hardware components such as RF dongles, specifically with the Realtek RTL2832U chipset, and a Raspberry Pi for the implementation of an SDR capable of real-time reception. The primary focus is on capturing FM mono/stereo audio and digital data using the Radio Data System (RDS) protocol.

## Features
- **RF Dongles and RTL2832U chipset:** The project leverages RF dongles equipped with the Realtek RTL2832U chipset for efficient radio frequency signal processing.

- **Raspberry Pi Integration:** A Raspberry Pi serves as the central processing unit, providing the computational power required for real-time radio signal analysis.

- **SDR Implementation:** The system is designed to function as a Software Defined Radio, allowing for the reception of FM mono/stereo audio and digital data.

- **RDS Protocol:** The project incorporates the RDS protocol for extracting additional information from the radio signals, enhancing the overall user experience.

- **Specialized Filtering Techniques:** To achieve optimal signal processing, the project applies specialized filtering techniques for the extraction, decimation, and demodulation of specific audio channels.

## Dependencies
- **RF Dongles with RTL2832U chipset**
- **Raspberry Pi**
- **C++ Compiler**
- **Python Interpreter**
- **SDR Libraries (e.g., GNU Radio)**
- **Additional Python Libraries (as specified in requirements.txt)**

## Usage
1. **Hardware Setup:** Connect the RF dongles with the RTL2832U chipset to the Raspberry Pi.

2. **Software Installation:**
   - Install necessary dependencies using the provided requirements.txt file.
   - Compile and execute the C++ components.

3. **Run the SDR Application:** Launch the Python script to initialize the real-time Software Defined Radio.

4. **Adjust Settings:** Modify parameters such as frequency, modulation type, and filtering options as needed.

5. **Explore Results:** Experience real-time reception of FM mono/stereo audio and digital data using the implemented SDR, with additional information decoded through the RDS protocol.

For detailed information and configuration options, please refer to the project documentation.

**Note:** Ensure proper permissions and configurations for the RF dongles and Raspberry Pi to enable seamless functionality of the Software Defined Radio.

## Acknowledgments

This project integrates C++ for microcontroller functionality and Python for data processing to create a Real-Time Software Defined Radio. Utilizing RF dongles, Realtek RTL2832U chipset, and a Raspberry Pi, the system achieves real-time reception of FM audio and digital data with RDS protocol. Specialized filtering techniques are applied for efficient signal processing. Check the project [report](doc/COE-3DY4%20Project%20Report.pdf) for specific components and constraints.
