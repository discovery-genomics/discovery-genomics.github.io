# Equipment

### Illumina NovaSeq 6000
- High throughput sequencer designed for speed and flexibility.
- Input: 48 samples in a 2 x 150 bp S4 flow cell
- Output after 44 hours: 2400-3000 GB
- The sequencer has an on-board hard drive. The hard drive is only for temporary storage before the data is transferred to a buffer server or BaseSpace. Having data saved in the hard drive that is different from the current run can decrease performance.

__System Software__
- Control computer operating system seems to be Windows.
- NovaSeq Control Software (NVCS): controls instrument operation, guides user through a run, displays statistics about the run.
- Real-Time Analysis (RTA): performs image analysis and base caling during run.
- Universal Copy Service (UCS): copies the output files from NVCS and RTA to output folder throughout a run. Also has the option to transfer to BaseSpace. If interrupted, UCS will attempt to reconnect multiple times and resume transfer of data.

__Storage__
- __Not sure how much onboard storage as of now__
- Compute Engine (CE) has storage and control computer has storage.

__Networking__
- Can configure output folder to be on a different network location.
- Illumina recommends 1 Gb connection between the sequencer and the buffer server.
- 200 MB/s/instrument for internal network uploads.
- 200 Mb/s/instrument for BaseSpace Sequence Hub uploads.

__Other Notes__
- Illumina Proactive monitoring services is on by default.

[Documentation](https://support.illumina.com/sequencing/sequencing_instruments/novaseq-6000/documentation.html)

### Illumina DRAGEN
