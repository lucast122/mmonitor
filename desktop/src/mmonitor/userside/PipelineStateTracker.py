import os
import json
from enum import Enum, auto

class PipelineStep(Enum):
    FLYE_ASSEMBLY = auto()
    MEDAKA_POLISH = auto()
    METABAT_BINNING = auto()
    CHECKM2_QUALITY = auto()
    BAKTA_ANNOTATION = auto()
    GTDBTK_TAXONOMY = auto()
    UPLOAD_TO_SERVER = auto()

class PipelineStateTracker:
    def __init__(self, output_dir, sample_name):
        self.output_dir = output_dir
        self.sample_name = sample_name
        
        # Create state file path in the sample directory
        sample_dir = os.path.join(output_dir, sample_name)
        self.state_file = os.path.join(sample_dir, "pipeline_state.json")
        
        # Ensure sample directory exists
        os.makedirs(sample_dir, exist_ok=True)
        
        self.state = self._load_state()

    @property
    def current_step(self):
        """Get the current pipeline step"""
        return self.state.get("current_step")

    def _load_state(self):
        """Load pipeline state from file"""
        if os.path.exists(self.state_file):
            try:
                with open(self.state_file, 'r') as f:
                    return json.load(f)
            except Exception as e:
                print(f"Error loading pipeline state: {e}")
        return {
            "steps_completed": [],
            "current_step": None,
            "step_outputs": {},
            "failed_step": None,
            "error_message": None,
            "pipeline_status": "not_started"  # Possible values: not_started, running, completed, failed
        }

    def _save_state(self):
        """Save pipeline state to file"""
        with open(self.state_file, 'w') as f:
            json.dump(self.state, f, indent=2)

    def initialize_pipeline_state(self):
        """Initialize pipeline state at the start of a run"""
        self.state["pipeline_status"] = "running"
        self.state["failed_step"] = None
        self.state["error_message"] = None
        self._save_state()

    def finalize_pipeline_state(self):
        """Finalize pipeline state at the end of a run"""
        if self.state.get("failed_step"):
            self.state["pipeline_status"] = "failed"
        else:
            self.state["pipeline_status"] = "completed"
        self._save_state()

    def _get_step_name(self, step) -> str:
        """Convert step to string name, handling both string and enum inputs"""
        if isinstance(step, PipelineStep):
            return step.name
        return str(step)

    def step_completed(self, step) -> bool:
        """Check if a step has been completed
        
        Args:
            step: Step name (string) or PipelineStep enum
        """
        step_name = self._get_step_name(step)
        return step_name in self.state["steps_completed"]
    
    # Alias for backward compatibility
    def is_step_complete(self, step) -> bool:
        """Alias for step_completed to maintain compatibility"""
        return self.step_completed(step)

    def mark_step_complete(self, step, output_files=None):
        """Mark a step as complete and save its output files
        
        Args:
            step: Step name (string) or PipelineStep enum
            output_files: Optional output files to store
        """
        step_name = self._get_step_name(step)
        if step_name not in self.state["steps_completed"]:
            self.state["steps_completed"].append(step_name)
        
        if output_files is not None:
            if "step_outputs" not in self.state:
                self.state["step_outputs"] = {}
            self.state["step_outputs"][step_name] = output_files
        
        self.state["current_step"] = None
        self.state["failed_step"] = None
        self.state["error_message"] = None
        self._save_state()

    def mark_step_failed(self, step, error_message):
        """Mark a step as failed with an error message
        
        Args:
            step: Step name (string) or PipelineStep enum
            error_message: Error message to store
        """
        step_name = self._get_step_name(step)
        self.state["failed_step"] = step_name
        self.state["error_message"] = error_message
        self.state["pipeline_status"] = "failed"
        self._save_state()

    def set_current_step(self, step: PipelineStep):
        """Set the current running step"""
        if step:
            self.state["current_step"] = step.name
            self.state["failed_step"] = None
            self.state["error_message"] = None
            self._save_state()

    def get_step_output(self, step):
        """Get the output files for a completed step
        
        Args:
            step: Step name (string) or PipelineStep enum
            
        Returns:
            The output path for the step, or None if step not completed
        """
        step_name = self._get_step_name(step)
        if not self.step_completed(step_name):
            return None
            
        output = self.state.get("step_outputs", {}).get(step_name)
        # If output is a dict, assume it's a legacy format and return the first value
        if isinstance(output, dict):
            return next(iter(output.values()), None)
        return output

    def set_step_output(self, step: PipelineStep, output_path: str):
        """Set the output path for a step

        Args:
            step (PipelineStep): The pipeline step
            output_path (str): Path to the step's output
        """
        if "step_outputs" not in self.state:
            self.state["step_outputs"] = {}
        self.state["step_outputs"][step.name] = output_path
        self._save_state()

    def get_state(self) -> dict:
        """Get the current pipeline state"""
        return self.state

    def get_pipeline_status(self):
        """Get the current status of the pipeline
        
        Returns:
            str: One of 'not_started', 'running', 'completed', or 'failed'
        """
        return self.state.get("pipeline_status", "not_started")

    def get_error_info(self):
        """Get information about any pipeline error
        
        Returns:
            tuple: (failed_step, error_message) or (None, None) if no error
        """
        return (self.state.get("failed_step"), self.state.get("error_message"))

    def reset(self):
        """Reset the pipeline state"""
        self.state = {
            "steps_completed": [],
            "current_step": None,
            "step_outputs": {},
            "failed_step": None,
            "error_message": None,
            "pipeline_status": "not_started"
        }
        self._save_state()
