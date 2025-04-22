# Structural Data Engineering of Uranium Species for New Metal Capture Agents

## Overview

This project investigates data-driven methods to analyze and design uranium-binding crown ether compounds, aiming to accelerate discovery for nuclear fuel reprocessing applications. Traditional trial-and-error approaches are augmented with a computational graph theory-based pipeline that parses crystallographic data to extract structural insights predictive of chelation efficacy.

## Background

Nuclear energy is a clean energy source, but effective uranium recovery from spent fuel is a major unsolved problem. Organic ligands such as crown ethers can bind uranium in its +VI oxidation state, but current ligand design is slow and inefficient.

This project leverages a **data engineering approach** to identify structural trends and predict effective uranium ligands—especially focusing on **macrocyclic crown ethers** and their interactions with the **uranyl ion (UO₂²⁺)**.

## Objectives

1. **Develop tools** to represent crystal structures as undirected graphs
2. **Identify and analyze** crown ether moieties in existing structural data
3. **Quantify geometric features** and correlate them with uranium-binding ability
4. **Expand the dataset** to explore a broader range of ligands with various donor atoms (O, N, S)
5. **Optimize data handling** using PostgreSQL for scalable, iterative analysis

## Project Timeline

| Phase | Duration        | Milestone                                          |
|-------|------------------|----------------------------------------------------|
| 1     | Month 1–3         | Research crown ether structures, develop pipeline |
| 2     | Month 4–6         | Apply vector-space and graph theory tools to 15-crown-5 data |
| 3     | Month 7–9         | Analyze trends, expand to more crown types        |
| 4     | Month 10–12       | Full dataset storage, visualization, and reporting |

## Methodology

- **Graph Theory**: Use BFS and DFS to convert molecules into undirected graphs where atoms are nodes and bonds are edges
- **3D Vector Space**: Analyze each molecule’s geometry using atomic positions
- **Feature Extraction**:
  - `ω_crown`: Planarity deviation
  - `Ψ_M`: Metal displacement from ring center
  - M–X bond lengths and coordination number
- **Storage & Querying**: Use PostgreSQL and other graph based databases to handle thousands of crystallographic records from the Cambridge Structural Database

## Deliverables

- Python scripts for parsing `.cif` files and generating graphs
- 3D visualizations of CE geometries
- Quantified database of CE structural features
- Hypothesis-tested trends on uranium chelation
- Public-facing website and interactive explorer:
	- Sample Website from previous project: [TFSI Research Site](https://tfsi-research.vercel.app/)
## References

1. [Pedersen, C.J., JACS 1967](https://pubs.acs.org/doi/abs/10.1021/ja01002a035)
2. [Golwankar et al., JACS 2024](https://pubs.acs.org/doi/10.1021/jacs.3c12075)
3. [Cambridge Structural Database](https://doi.org/10.1107/S0108768102003890)
4. [Depth-First Search – GeeksForGeeks](https://www.geeksforgeeks.org/depth-first-search-or-dfs-for-a-graph/)
5. [PostgreSQL](https://www.postgresql.org/)

## Acknowledgments

This project is being pursued under the mentorship of James Blakemore and his group.

---

> _Author: Tej Gumaste_  
> _Major: Computer Science | Minor: Mathematics_  
> _University of Kansas_

