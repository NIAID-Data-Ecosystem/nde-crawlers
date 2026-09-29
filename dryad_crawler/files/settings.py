LOG_LEVEL = "INFO"
DOWNLOAD_DELAY = 0.5
AUTOTHROTTLE_ENABLED = True
AUTOTHROTTLE_TARGET_CONCURRENCY = 1.0
AUTOTHROTTLE_DEBUG = True
ROBOTSTXT_OBEY = True
# Dryad allows anonymous clients 30 API requests per minute per IP; a fixed 2.5 s delay keeps us at 24
DOWNLOAD_SLOTS = {"dryad-api": {"concurrency": 1, "delay": 2.5, "jitter": 0}}
HTTPCACHE_ENABLED = True
HTTPCACHE_EXPIRATION_SECS = 0
HTTPCACHE_DIR = "/cache"
HTTPCACHE_POLICY = "middlewares.CachePolicy"
