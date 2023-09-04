from dataclasses import dataclass

@dataclass
class Diagram:
    paired: set[tuple[int, int]]
    unpaired: set[int]

@dataclass
class DiagramEnsemble:
    ker: Diagram
    cok: Diagram
    dom: Diagram
    cod: Diagram 
    im: Diagram 
    rel: Diagram