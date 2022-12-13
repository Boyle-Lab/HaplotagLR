import nox


@nox.session
def tests(session):
    """
    Args:
        session:
    """
    session.install('pytest')
    session.run('pytest')


@nox.session
def lint(session):
    """
    Args:
        session:
    """
    session.install('flake8')
    session.run('flake8', '--import-order-style', 'google')
