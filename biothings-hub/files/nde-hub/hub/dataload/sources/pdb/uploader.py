from hub.dataload.nde import NDESourceUploader

# PDB is only in scope for NIAID: keep entries with an infectious agent, or
# funded by NIAID (ROR 043z4tv69).
NIAID_ROR = "https://ror.org/043z4tv69"


class PDB_Uploader(NDESourceUploader):
    name = "pdb"

    def post_process(self, doc):
        in_scope = doc.get("infectiousAgent") or any(
            funder.get("identifier") == NIAID_ROR
            for funding in doc.get("funding", [])
            for funder in funding.get("funder", [])
            if isinstance(funder, dict)
        )
        return doc if in_scope else None
