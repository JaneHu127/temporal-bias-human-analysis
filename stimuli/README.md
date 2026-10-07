# Stimuli

This folder contains the stimulus materials used in the temporal reconstruction experiment, including the original video stimuli and the image triplets used in the temporal-location reconstruction task.

## Folder structure

```text
stimuli/
├── videos/
└── triplets/
```

- `videos/` contains the video clips used in the experiment.
- `triplets/` contains the image stimuli used in the temporal reconstruction trials.

## Video stimuli

The `videos/` folder contains the video clips used in the experiment.

Frames extracted from these videos were used to construct the triplet stimuli for the temporal reconstruction task. Each triplet consisted of a **start frame**, a **target frame**, and an **end frame**.

## Triplet stimuli

The triplet stimuli were organized into three interval-width conditions:

- **12 s**
- **60 s**
- **300 s**

For each interval-width condition, **25 pairs of start and end frames** were selected.

The stimulus sets are indexed from **1 to 75**:

- **1–25**: 12 s interval-width condition
- **26–50**: 60 s interval-width condition
- **51–75**: 300 s interval-width condition

Each numbered stimulus set corresponds to one pair of start and end frames.

## Target temporal locations

For each start–end frame pair, four target frames were selected at four relative temporal locations within the interval:

| Target condition | Relative temporal location |
| --- | --- |
| T20 | 20% |
| T40 | 40% |
| T60 | 60% |
| T80 | 80% |

The four target frames were used in separate trials.

Thus, each start–end frame pair generated four triplets:

- Start + Target 20% + End
- Start + Target 40% + End
- Start + Target 60% + End
- Start + Target 80% + End

For each interval-width condition, there were:

- 25 start–end frame pairs
- 4 target locations per pair
- 100 triplets

Across all three interval-width conditions, this resulted in **75 start–end frame pairs** and **300 triplets** in total.

## Trial structure

In each temporal reconstruction trial, three images were presented simultaneously:

1. **Start frame** – the beginning of the corresponding temporal interval.
2. **Target frame** – a frame located at 20%, 40%, 60%, or 80% of the interval.
3. **End frame** – the end of the corresponding temporal interval.

Participants were asked to estimate the temporal position of the target frame relative to the start and end frames.

The interval defined by the start and end frames was normalized to a **0–100 scale**, where:

- **Start frame = 0**
- **End frame = 100**

The target frame therefore corresponded to position 20, 40, 60, or 80 on this normalized temporal axis.

## Organization of the triplet images

For convenience, the materials associated with each start–end frame pair are displayed together as a six-image set:

**Start | Target 20% | Target 40% | Target 60% | Target 80% | End**

The leftmost and rightmost images are the shared **start** and **end** frames, respectively. The four images in between are the four possible **target frames** associated with that cue pair.

Importantly, participants did **not** see all six images in a single trial. Each experimental trial contained only one start frame, one target frame, and one end frame.
