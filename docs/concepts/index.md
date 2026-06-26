---
title: Concepts
description: High-level genome assembly concepts that explain why ChromoSort workflows work.
---

# Concepts

These pages are the educational frame for ChromoSort. They teach the ideas a
new student or researcher should understand before choosing a command: what the
evidence means, why the action is reasonable, and where the action can go
wrong.

| Concept | Use it when |
| --- | --- |
| [How to interpret dot plots]({{ '/dot-plots/' | relative_url }}) | You need the visual grammar for reference/query alignment plots before deciding what action is justified. |
| [How to decide when to fix a contig]({{ '/concepts/fixing-contigs/' | relative_url }}) | You need to decide whether a suspicious alignment pattern is strong enough to justify splitting or cutting sequence. |
| [How scaffolding works]({{ '/concepts/scaffolding/' | relative_url }}) | You need to understand how ordered contigs become chromosome-scale scaffold records and why Ns are placeholders, not discovered sequence. |
| [How gap filling works]({{ '/concepts/gap-filling/' | relative_url }}) | You need to understand when graph or read evidence can replace scaffold Ns with sequence, and why ambiguous fills should stay unresolved. |

After the concept is clear, move to the [Guides]({{ '/guides/' | relative_url }})
for ChromoSort-specific inputs, commands, reports, and reviewed application
steps.
