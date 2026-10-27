import numpy as np
from astropy import units as u

from tardis.io.logger.logger import TARDISLogger
from tardis.util.environment import Environment
from tardis.workflows.workflow_logger import WorkflowLogger


def test_log_plasma_state_with_specific_log_level(monkeypatch, caplog):
    monkeypatch.setattr(
        Environment, "allows_widget_display", staticmethod(lambda: True)
    )
    tardis_logger = TARDISLogger()
    tardis_logger.configure_logging(
        "INFO", {"debug": {}}, specific_log_level=True
    )
    workflow_logger = WorkflowLogger.__new__(WorkflowLogger)

    try:
        workflow_logger.log_plasma_state(
            t_rad=np.array([1e4, 9e3]) * u.K,
            dilution_factor=np.array([0.5, 0.4]),
            t_inner=1e4 * u.K,
            next_t_rad=np.array([1.1e4, 9.5e3]) * u.K,
            next_dilution_factor=np.array([0.55, 0.45]),
            next_t_inner=1.05e4 * u.K,
        )
    finally:
        tardis_logger.configure_logging("INFO", {}, specific_log_level=False)

    assert "Current t_inner" in caplog.records[-1].getMessage()
