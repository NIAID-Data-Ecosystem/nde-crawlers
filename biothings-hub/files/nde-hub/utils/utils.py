import functools
import time
import traceback
from config import logger

def retry(retry_num, retry_sleep_sec):
    """
    retry help decorator.
    :param retry_num: the retry num; retry sleep sec
    :return: decorator
    """
    def decorator(func):
        """decorator"""
        # preserve information about the original function, or the func name will be "wrapper" not "func"
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            """wrapper"""
            for attempt in range(retry_num):
                try:
                    return func(*args, **kwargs)  # should return the raw function's return value
                except Exception as err:   
                    logger.error(err)
                    logger.error(traceback.format_exc())
                    time.sleep(retry_sleep_sec)
                logger.info("Retrying failed func %s. Trying attempt %s of %s.", func, attempt + 1, retry_num)
            logger.error("func %s retry failed", func)
            raise Exception('Exceed max retry num: {} failed'.format(retry_num))

        return wrapper

    return decorator

def check_schema(record):
    """
    check record before inserting into MongoDb
    :param record: record to be checked
    """
    
    assert isinstance(record, dict), "record is not a dict"
    assert record.get("_id"), "_id is None"
    assert record.get("@type"), "@type is None"
    assert record.get("includedInDataCatalog"), "includedInDataCatalog is None"
    if coa := record.get("conditionsOfAccess"):
        enum = ["Open", "Restricted", "Closed", "Embargoed"]
        assert coa in enum, "%s is not a valid conditionsOfAccess. Allowed conditionsOfAccess: %s" % (coa, enum)