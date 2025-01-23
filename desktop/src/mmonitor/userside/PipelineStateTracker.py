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
        self.state_file = os.path.join(output_dir, sample_name, "pipeline_state.json")
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
            "error_message": None
        }

    def _save_state(self):
        """Save pipeline state to file"""
        os.makedirs(os.path.dirname(self.state_file), exist_ok=True)
        with open(self.state_file, 'w') as f:
            json.dump(self.state, f, indent=2)

    def step_completed(self, step: PipelineStep) -> bool:
        """Check if a step has been completed"""
        return step.name in self.state["steps_completed"]

    def mark_step_complete(self, step: PipelineStep, output_files=None):
        """Mark a step as complete and save its output files"""
        if output_files:
            self.state["step_outputs"][step.name] = output_files
        if step.name not in self.state["steps_completed"]:
            self.state["steps_completed"].append(step.name)
        self.state["current_step"] = None
        self.state["failed_step"] = None
        self.state["error_message"] = None
        self._save_state()

    def mark_step_failed(self, step: PipelineStep, error_message: str):
        """Mark a step as failed with an error message"""
        if step:
            self.state["failed_step"] = step.name
            self.state["error_message"] = error_message
            self.state["current_step"] = None
            self._save_state()

    def set_current_step(self, step: PipelineStep):
        """Set the current running step"""
        if step:
            self.state["current_step"] = step.name
            self.state["failed_step"] = None
            self.state["error_message"] = None
            self._save_state()

    def get_step_output(self, step: PipelineStep) -> str:
        """Get the output path for a completed step

        Args:
            step (PipelineStep): The pipeline step

        Returns:
            str: The output path for the step, or None if step not completed
        """
        if not self.step_completed(step):
            return None
            
        output = self.state.get("step_outputs", {}).get(step.name)
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

    def reset(self):
        """Reset the pipeline state"""
        self.state = {
            "steps_completed": [],
            "current_step": None,
            "step_outputs": {},
            "failed_step": None,
            "error_message": None
        }
        self._save_state()
