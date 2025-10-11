
#ifndef SOLVER_METADATA_EXTRACT_H
#define SOLVER_METADATA_EXTRACT_H

void solver_metadata_extract(struct solver_metadata *metadata, const Quadrature* quad, const Mesh* mesh,
	                         const XSection* xs, const SourceParams* srcPar, int pn=0);

#endif
