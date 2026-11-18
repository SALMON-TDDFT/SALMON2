"""Session-local native factory ownership; no loading at import."""
from .errors import RuntimeFailure

class _SessionGroup:
    """Session-local model owner, including any returned but unadmitted states."""
    def __init__(self,model):
        self.model=model
        self.states=[]
    def new_state(self,point_id):
        state=self.model.new_state(point_id)
        self.states.append(state)
        return state
    def close(self):
        errors=[]
        for state in self.states:
            try:state.close()
            except Exception as exc:errors.append(str(exc))
        try:self.model.close()
        except Exception as exc:errors.append(str(exc))
        if errors:raise RuntimeFailure('shared native group cleanup failed: '+str(errors))
