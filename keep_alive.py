"""
Selenium-based keep-alive for Streamlit Community Cloud.

Streamlit Cloud hibernates apps after ~12h of inactivity. A simple HTTP
ping (curl) only reaches the reverse proxy — it never establishes the
WebSocket session that Streamlit counts as "real" user activity.

This script uses a headless Chrome browser to:
  1. Navigate to the app URL (establishing a full browser session).
  2. Detect and click the "Yes, get this app back up!" button if the
     app is already sleeping.
  3. Hold the session open long enough for the WebSocket connection to
     register as genuine activity, resetting the inactivity timer.

Intended to be run by a GitHub Actions cron job every 6 hours.
"""

import os
import sys
import time

from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import WebDriverWait
from webdriver_manager.chrome import ChromeDriverManager


def main():
    url = os.environ.get("APP_URL")
    if not url:
        print("::error::APP_URL environment variable is not set.")
        print("Add it as a repository secret: Settings → Secrets → Actions → APP_URL")
        sys.exit(1)

    print(f"🔔 Pinging: {url}")
    print(f"⏰ Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}")

    # ── Headless Chrome options (CI-friendly) ────────────────────────
    opts = Options()
    opts.add_argument("--headless=new")
    opts.add_argument("--no-sandbox")
    opts.add_argument("--disable-dev-shm-usage")
    opts.add_argument("--disable-gpu")
    opts.add_argument("--window-size=1920,1080")

    driver = webdriver.Chrome(
        service=Service(ChromeDriverManager().install()),
        options=opts,
    )

    try:
        driver.get(url)

        # ── Phase 1: Check if app is sleeping ────────────────────────
        try:
            wake_btn = WebDriverWait(driver, 20).until(
                EC.element_to_be_clickable(
                    (By.XPATH, "//button[contains(., 'Yes, get this app back up')]")
                )
            )
            print("😴 App is sleeping — clicking wake-up button...")
            wake_btn.click()

            # Wait for the main Streamlit container to appear (app loaded)
            WebDriverWait(driver, 60).until(
                EC.presence_of_element_located(
                    (By.CSS_SELECTOR, "[data-testid='stAppViewContainer']")
                )
            )
            print("✅ App successfully woken up!")

        except Exception:
            # No wake-up button found within 20s → app is already awake
            print("☕ App is already awake — establishing session...")

        # ── Phase 2: Hold session to register as real activity ───────
        # Streamlit's WebSocket connects when the page loads.  Keeping
        # the browser open for a few seconds ensures the connection is
        # fully established and counted as "active" by the platform.
        time.sleep(15)
        print("✅ Session held for 15 seconds — activity registered.")
        print(f"📄 Final page title: {driver.title}")

    except Exception as e:
        print(f"::warning::Keep-alive encountered an issue: {e}")
        # Don't fail the workflow — next scheduled run will retry
        sys.exit(0)

    finally:
        driver.quit()


if __name__ == "__main__":
    main()
