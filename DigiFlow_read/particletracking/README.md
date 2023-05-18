## Particle Tracking Codes
Tools are provided here for the purposes of tracking surface particles and producing a coherent .mat structure describing them, either via post-processing of DigiFlow PTV output, or from digiflow images. See the following flow chart to identify the process:

```mermaid
flowchart TD
    A[Is Float motion captured from above?] -->|Yes| C[DigiFlow PTV]


